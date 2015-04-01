#include "Viewer.h"
#include "ui_Viewer.h"
#include <QPainter>
#include <QDebug>
#include <QTimer>
#include <QKeyEvent>
#include <QMouseEvent>
#include <QFileDialog>

#include "libfastmarching.h"

Viewer::Viewer(QWidget *parent) : QGLWidget(parent), ui(new Ui::Viewer)
{
	QGLFormat glf = QGLFormat::defaultFormat();
	glf.setSampleBuffers(true);
	this->setFormat(glf);

    ui->setupUi(this);
	setFocusPolicy(Qt::ClickFocus);
}

Viewer::~Viewer()
{
    delete ui;
}

void Viewer::paintEvent(QPaintEvent *)
{
    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing);
    painter.fillRect(rect(),Qt::white);

    for(int i = 0; i < images.size(); i++)
    {
        MyQImage & image = images[i];
        QPoint pos(image.x, image.y);

        painter.translate(pos);
        painter.setOpacity(image.opacity);
        painter.drawImage(0,0, image);
        painter.setOpacity(1.0);

        //MyPolygon & polygon = polys[i];
        //polygon.draw( painter, pos + QPoint(0, this->height()), image.width(), image.height() );

		painter.translate(-pos);
    }

    painter.save();
    painter.translate(2,2);
    painter.setPen(QPen(QColor(0,0,0,128), 20, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
    painter.drawPoints(points);
    painter.restore();
    painter.setPen(QPen(Qt::blue, 20, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
    painter.drawPoints(points);

    for(auto poly: polys)
    {
        painter.drawPolygon(poly);
    }
}

template<typename QPolygonType>
QPolygonType resamplePolygon(QPolygonType points, int count = 100){
	QPainterPath path;
	path.addPolygon(points);
	auto pathLen = path.length();
	auto stepSize = pathLen / count;
    QPolygonType newPoints;
	for(int i = 0; i < count; i++)
    {
        QPolygon::value_type p = path.pointAtPercent(path.percentAtLength( stepSize * i )).toPoint();
        newPoints << p;
    }
	return newPoints;
}

template<typename QPolygonType>
QPolygonType smoothPolygon( QPolygonType points, int iterations = 1 ){
	for(int i = 0; i < iterations; i++)
	{
        QPolygonType newPoints;

        for(int p = 0; p < points.size(); p++)
        {
			int s = p-1, t = p+1;
            if(s < 0) s = 0;
            if(t > points.size()-1) t = points.size()-1;
            auto prev = points[s];
            auto next = points[t];
			newPoints << (prev + next) * 0.5;
		}

		points = newPoints;
	}
	return points;
}

void Viewer::mousePressEvent(QMouseEvent *event)
{
    if(event->buttons().testFlag(Qt::LeftButton))
    {
        if(event->modifiers().testFlag(Qt::ControlModifier))
        {
            points << event->pos();
        }

        if(event->modifiers().testFlag(Qt::ShiftModifier))
        {
            points.clear();
        }
    }

    update();
}

void Viewer::keyPressEvent(QKeyEvent *event)
{
    if(event->key() == Qt::Key_I)
    {
        images.clear();
        images << MyQImage(QImage(QFileDialog::getOpenFileName(0,"Load image", "", "*.png")));

        points.clear();
    }

    if(event->key() == Qt::Key_S)
    {
        if(points.size() < 2) return;

        auto shape = images.front();
        auto start_point = points.at(points.size()-2).toPoint();
        auto end_point = points.back().toPoint();

        libfastmarching::Path fm2path;
        libfastmarching::fm2star(shape, start_point, end_point, fm2path);

        // Convert
        QPolygon path;
        for(auto p : fm2path) path << QPoint(p[0], p[1]);

        // Post-process path
        path = smoothPolygon(resamplePolygon(path,shape.width() * 0.25), 2);

        // Visualize
        {
            QImage strokeImage(shape.width(), shape.height(), QImage::Format_ARGB32_Premultiplied);
            strokeImage.fill(Qt::transparent);
            QPainter painter(&strokeImage);
            painter.setRenderHint(QPainter::Antialiasing);

            // Stroke color
            QColor penColor;
            penColor.setHslF(double(rand()) / RAND_MAX * 0.25, 1, 0.5);
            painter.setPen(QPen(penColor, 40, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));

            // Draw stroke
            painter.drawPolyline(path);

            images << MyQImage(strokeImage);
        }
    }

    update();
}
