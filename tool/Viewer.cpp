#include "Viewer.h"
#include "ui_Viewer.h"
#include <QPainter>
#include <QDebug>
#include <QTimer>
#include <QKeyEvent>

#include "MarchingSquares.h"

Viewer::Viewer(QWidget *parent) : QGLWidget(parent), ui(new Ui::Viewer)
{
	QGLFormat glf = QGLFormat::defaultFormat();
	glf.setSampleBuffers(true);
	this->setFormat(glf);

    ui->setupUi(this);
	setFocusPolicy(Qt::ClickFocus);

    // Experiment
    demo();
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

    for(int i = 0; i < polys.size(); i++)
    {
        MyQImage & image = images[i];
        QPoint pos(image.x, image.y);

        painter.translate(pos);
        painter.setOpacity(image.opacity);
        painter.drawImage(0,0, image);
        painter.setOpacity(1.0);

        MyPolygon & polygon = polys[i];
		polygon.draw( painter, pos + QPoint(0, this->height()), image.width(), image.height() );

		painter.translate(-pos);
    }
}

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

void Viewer::demo()
{
    // Draw some letters
    int x = 0;
    QString word = "monkey";
    for(auto letter : word)
    {
        QImage img(256, 256, QImage::Format_ARGB32_Premultiplied);
        img.fill(qRgba(0,0,0,0));

        // Draw letter
        {
            QPainter painter( &img );
            painter.setRenderHint(QPainter::Antialiasing, false);
            //painter.drawRect(img.rect());
            painter.setFont(QFont("Arial", 200, QFont::Bold));
            painter.drawText( img.rect(), Qt::AlignCenter,  QString("%1").arg(letter));
        }

        // Convert to double matrix
        Eigen::MatrixXd img_matrix = Eigen::MatrixXd::Zero(img.height(), img.width());
        for(int y = 0; y < img.height(); y++){
            for(int x = 0; x < img.width(); x++){
                if(img.pixel(x,y) > 0)
                {
                    img_matrix(y,x) = 1.0;
                }
            }
        }

        // Extract outermost boundary
        QPolygonF poly;
        for(auto pixel : MarchingSquares(img_matrix, 1.0).march())
            poly << QPointF( pixel.x(), pixel.y() );


		poly = resamplePolygon(poly, 100);
		//poly = smoothPolygon(poly, 2);

        polys.push_back( MyPolygon( poly, QString("letter_%1").arg(letter) ) );

        // Display extracted boundary
        {
            QPainter painter(&img);
            painter.setPen(QPen(Qt::blue, 3));
            painter.drawPolygon(poly);
        }

        images.push_back( MyQImage(img, x, 0, 0.1) );

        x += img.width();
    }
}

void Viewer::keyPressEvent(QKeyEvent *event)
{
	if(event->key() == Qt::Key_Space){
		QTimer * timer = new QTimer;
		connect(timer, &QTimer::timeout, [=](){
			for(auto & poly : polys) poly.deform();
			update();
		});
		timer->start( 20 );
	}
}
