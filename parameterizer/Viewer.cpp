#include "Viewer.h"
#include "ui_Viewer.h"
#include <QPainter>
#include <QKeyEvent>
#include <QFileInfo>
#include <QTextStream>
#include <iostream>

using namespace std;

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

QPolygonF smoothPolygon( QPolygonF points, int iterations = 1 ){
    for(int i = 0; i < iterations; i++)
    {
        QPolygonF newPoints;
        points.removeLast();

        for(int p = 0; p < points.size(); p++){
            int s = p-1, t = p+1;
            if(s < 0) s = points.size() - 1;
            if(t > points.size()-1) t = 0;
            QPointF prev = points[s];
            QPointF next = points[t];
            newPoints << (prev + next) * 0.5;
        }

        newPoints << newPoints.front();
        points = newPoints;
    }
    return points;
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
        contours.push_back(contour);
    }
}

Viewer::~Viewer()
{
	delete ui;
}

void Viewer::paintEvent(QPaintEvent *)
{
    QPainter painter(this);

    for (auto c : contours)
    {
        painter.drawPath(c);
    }

    //painter.fillRect(rect(), QColor(255, 0, 0));

    //painter.setPen(QPen(Qt::black, 2));
    //painter.drawText(QPoint(rect().center()), message);

	painter.setPen(QPen(Qt::green, 4));
    for (auto p : paths){
		painter.drawPath(p);
    }
}
void Viewer::keyPressEvent(QKeyEvent *event)
{
    /*if (event->key() == Qt::Key_Space)
	{
		message = QString("Hello! a number %1 and %2").arg(3.14).arg(" pi");
		update();
    }*/
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

		QPainterPath path;
		path.moveTo(event->pos());
		paths.push_back(path);
	}
}
void Viewer::mouseReleaseEvent(QMouseEvent *event)
{
	if (event->button() == Qt::LeftButton && scribbling) {
		scribbling = false;

        if(paths.size())
        {
            auto & path = paths.back();

            QPainterPath smoothPath;
            QPolygonF poly = path.toFillPolygons().front();

            if(!poly.empty())
            {
                poly = resamplePolygon(poly, 10);

                smoothPath.moveTo(poly.first());

                for(int i = 1; i < (int)poly.size() + 1; i += 3)
                {
                    smoothPath.quadTo(poly[i], poly[i+1]);
                }
            }
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

        p.lineTo(QPoint(x, y));

		update();
	}
}

void Viewer::reduceDimention()
{
    qreal test;
    for (auto p : paths)
    {
        test = p.angleAtPercent(0.5);
        cout<<"test:"<<test<<endl;
    }
}
