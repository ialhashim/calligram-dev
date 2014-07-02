#include <math.h>
#define _USE_MATH_DEFINES

#include "Viewer.h"
#include "ui_Viewer.h"
#include <QPainter>
#include <QKeyEvent>
#include <QFileInfo>
#include <QTextStream>
#include <iostream>

using namespace std;

// Eigen library
#include <Eigen/Core>
#include <Eigen/Geometry>
using namespace Eigen;
typedef Vector2d Vector2;

// Intersection library
#include "intersections.h"

// Helper functions:
QPolygonF resamplePolygon(QPolygonF points, int count = 100){
	QPainterPath path;
	path.addPolygon(points);
	auto pathLen = path.length();
	auto stepSize = pathLen / count;
	QPolygonF newPoints;
	for(int i = 0; i < count; i++)
		newPoints << path.pointAtPercent(path.percentAtLength( stepSize * i ));
	newPoints << path.pointAtPercent(0);
	return newPoints;
}
QPolygonF smoothPolygon( QPolygonF points, int iterations){
	for(int i = 0; i < iterations; i++)
	{
		QPolygonF newPoints;

		for(int p = 0; p < points.size(); p++)
		{
			int s = p-1, t = p+1;
			if(s < 0) s = 0;
			if(t > points.size()-1) t = points.size() - 1;
			QPointF prev = points[s];
			QPointF next = points[t];
			newPoints << (prev + next) * 0.5;
		}

		points = newPoints;
	}
	return points;
}

bool isIntersect( const QLineF & line1, const QPainterPath & path, QPointF & isect )
{
	QPolygonF polygon = path.toFillPolygon();

	QMap<double, QPointF> isections;

	for(size_t i = 0; i < polygon.size(); i++)
	{
		int j = i + 1;
		if(i == polygon.size()-1) j = 0;

		QLineF line2( polygon[i], polygon[j] );

		//intersection2D::Point isect1, isect2;
		//intersection2D::Segment s1();
		//intersection2D::Segment s2();
		//int isIntersect = intersection2D::intersect2D_2Segments( s1, s2, &isect1, &isect2 );

		float x, y;

		int isIntersect = get_line_intersection(line1.x1(), line1.y1(), line1.x2(), line1.y2(), 
												line2.x1(), line2.y1(), line2.x2(), line2.y2(), &x, &y);

		if( isIntersect > 0 )
		{
			QPointF cur_isect = QPointF( x, y );
			QPointF delta = (cur_isect - line1.p1());
			double distance = Vector2(delta.x(), delta.y()).norm();

			isections[ distance ] = cur_isect;
		}
	}

	if( isections.empty() )
		return false;
	else
	{
		isect =  isections[ isections.keys().front() ];
		return true;
	}
}

Viewer::Viewer(QWidget *parent) :
	QWidget(parent),
	ui(new Ui::Viewer)
{
	ui->setupUi(this);

	setFocusPolicy(Qt::ClickFocus);

	QString defaultFile = "../points.txt";

	QPolygonF points;

	// Load points from file
	if(QFileInfo(defaultFile).exists()){
		QFile file( defaultFile );
		file.open(QFile::ReadOnly | QFile::Text);
		QTextStream in(&file);
		QStringList lines = in.readAll().split('\n');
		for( auto line : lines ){
			auto pointCoord = line.split(",");
			if(pointCoord.size() < 2) continue;
			points << QPointF( pointCoord[0].toDouble(), pointCoord[1].toDouble() );
		}
		points << points.first();
	}
	else
		return;

	QPainterPath contour;
	bool isFirst = true;

	for(auto p : points)
	{
		float x = p.x();
		float y = p.y();

		if(isFirst) {
			contour.moveTo(QPoint(x,y));
			isFirst = false;
		}

		contour.lineTo(QPoint(x,y));
	}
	contours.push_back(contour);
}

Viewer::~Viewer()
{
	delete ui;
}

void Viewer::paintEvent(QPaintEvent *)
{
	QPainter painter(this);
	painter.setRenderHint(QPainter::Antialiasing);

	if(contours.empty()) return;

	for (auto & c : contours)
	{
		painter.drawPath(c);
	}

	painter.setPen(QPen(Qt::green, 4));
	for (auto & path : paths)
	{
		painter.drawPolyline(path);
	}

	for (auto & polygon : paths)
	{
		QVector<double> distancesTop, distancesBottom;

		QPainterPath path;
		path.addPolygon(polygon);

		int samplesCount = 60;

		auto pathLen = path.length();
		auto stepSize = pathLen / samplesCount;

		QPolygonF tangents;

		int pointSize = 4;
		int lineWidth = 2;

		int outside = 200;

		for(int i = 0; i <= samplesCount; i++)
		{
			QPointF pf, qf, mf;

			// Current point
			pf = path.pointAtPercent(path.percentAtLength( stepSize * i ));
			
			// Next point
			if(i < samplesCount-1) 
			{
				qf = path.pointAtPercent(path.percentAtLength( stepSize * (i+1) ));
			}
			else
			{
				mf = path.pointAtPercent(path.percentAtLength( stepSize * (i-1) ));
				qf = (1 * (pf-mf)) + pf;
			}

			painter.setPen(QPen(Qt::red, pointSize));
			painter.drawPoint(pf);

			Vector2 p(pf.x(),pf.y()), q(qf.x(),qf.y());

			if( (q-p).norm() < 1e-5 ) continue; // Segment too short

			Vector2 t( (q-p).normalized() );

			q = p + t * (p-q).norm();

			painter.setPen(QPen(Qt::blue, lineWidth * 0.5));
			painter.drawLine( p.x(), p.y(), q.x(), q.y() );

			Rotation2D<double> rot( M_PI_2 );
			Vector2 w = rot * t;
			Vector2 u1 = p + w * outside;
			Vector2 u2 = p + w * -outside;

			painter.setPen(QPen(QColor(255,255,0,80), lineWidth));
			painter.drawLine( p.x(), p.y(), u1.x(), u1.y());
			painter.drawLine( p.x(), p.y(), u2.x(), u2.y());

			QPainterPath cont = contours.front();
			QPointF isect;
			QLineF raySegment1( p.x(), p.y(), u1.x(), u1.y() );
			QLineF raySegment2( p.x(), p.y(), u2.x(), u2.y() );

			if( isIntersect(raySegment1, cont, isect) )
			{
				painter.setPen(QPen(Qt::red, 5));
				painter.drawPoint( isect );

				distancesBottom.push_back( (p - Vector2(isect.x(), isect.y())).norm() );
			}

			if( isIntersect(raySegment2, cont, isect) )
			{
				painter.setPen(QPen(Qt::green, 5));
				painter.drawPoint( isect );

				distancesTop.push_back( (p - Vector2(isect.x(), isect.y())).norm() );
			}
		}

		// Draw the function
		QPainterPath f1, f2;

		int width = 300;
		int height = 100;

		double minf1 = *std::min_element(distancesBottom.begin(), distancesBottom.end());
		double minf2 = *std::min_element(distancesTop.begin(), distancesTop.end());
		
		double maxf1 = *std::max_element(distancesBottom.begin(), distancesBottom.end());
		double maxf2 = *std::max_element(distancesTop.begin(), distancesTop.end());

		QLineF l1(50,50,50+width,50);
		QLineF l2 = l1.translated( QPointF(width,0) );

		f1.moveTo(l1.p1());
		f2.moveTo(l2.p1());

		size_t N = std::min(distancesTop.size(), distancesBottom.size());

		for(size_t i = 0; i < N; i++)
		{
			double t = double(i) / (distancesTop.size()-1);

			// Normalized
			//double v1 = (distancesTop[i] - minf1) / (maxf1-minf1);
			//double v2 = (distancesBottom[i] - minf2) / (maxf2-minf2);

			double v1 = distancesTop[i] / 50;
			double v2 = distancesBottom[i] / 50;

			QPointF p1 = l1.pointAt(t);
			QPointF p2 = l2.pointAt(t);

			f1.lineTo( p1 - QPoint(0,v1 * height)  );
			f2.lineTo( p2 - QPoint(0,v2 * height)  );
		}

		painter.translate(0, 600);

		painter.setPen(QPen(Qt::green, 2));
		painter.drawPath(f1);

		painter.setPen(QPen(Qt::red, 2));
		painter.drawPath(f2);
	}
}
void Viewer::keyPressEvent(QKeyEvent *event)
{
	if (event->key() == Qt::Key_Space)
	{
		reduceDimention();
	}
}

void Viewer::mousePressEvent(QMouseEvent *event)
{
	if (event->button() == Qt::LeftButton) {
		lastPoint = event->pos();
		scribbling = true;

		QPolygonF path;
		paths.push_back(path);
	}
	else
	{
		paths.clear();
		update();
	}
}
void Viewer::mouseReleaseEvent(QMouseEvent *event)
{
	if (event->button() == Qt::LeftButton && scribbling) {
		scribbling = false;

		if(!paths.empty())
		{
			auto & path = paths.back();
			path = smoothPolygon(path, 10);
		}
	}

	update();
}

void Viewer::mouseMoveEvent(QMouseEvent *event)
{
	if ((event->buttons() & Qt::LeftButton) && scribbling)
	{
		int prev_y = lastPoint.y();
		int prev_x = lastPoint.x();

		int y = event->pos().y();
		int x = event->pos().x();

		auto & p = paths.back();

		p << QPoint(x, y);

		update();
	}
}

void Viewer::reduceDimention()
{
	qreal test;
	for (auto p : paths)
	{
		//test = p.angleAtPercent(0.5);
		cout<<"test:"<<test<<endl;
	}
}
