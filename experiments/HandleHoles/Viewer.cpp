#include "Viewer.h"
#include "ui_Viewer.h"

#include <QPainter>
#include <QDebug>
#include <QTimer>
#include <QKeyEvent>
#include <QMouseEvent>
#include <QFileDialog>
#include <QTextStream>
#include <QString>
#include <QFile>
#include <QStringList>
#include <QLineF>
#include <QtAlgorithms>
#include <QVector3D>

#include "../ImageToVector/imagetovector.h"

Q_DECLARE_METATYPE(QPainterPath)

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

        painter.translate(-pos);
    }

    painter.setPen(QPen(Qt::red, 3));
    painter.setBrush(Qt::transparent);
    for(auto p : polys)
    {
        painter.drawPolygon(p);
    }

    int pointSize = 10;

    painter.save();
    painter.translate(0.3 * pointSize,0.3 * pointSize);
    painter.setPen(QPen(QColor(0,0,0,128), pointSize, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
    painter.drawPoints(points);
    painter.restore();
    painter.setPen(QPen(Qt::blue, pointSize, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
    painter.drawPoints(points);

    // DEBUG:
    for(QVariant v : debugs)
    {
        painter.setPen(QPen(Qt::blue, 2));
        if(QString(v.typeName()) == "QPolygonF") painter.drawPolygon(v.value<QPolygonF>());
        if(QString(v.typeName()) == "QLineF") painter.drawLine(v.value<QLineF>());
        if(QString(v.typeName()) == "QPointF") {painter.setPen(QPen(Qt::blue, 8)); painter.drawPoint(v.value<QPointF>());}
        if(QString(v.typeName()) == "QPainterPath") painter.drawPath(v.value<QPainterPath>());
    }
}

void Viewer::mousePressEvent(QMouseEvent *event)
{
    if(event->buttons().testFlag(Qt::LeftButton))
    {
        if(event->modifiers().testFlag(Qt::ControlModifier))
        {
            points << event->pos();

            // orthognal vector
            if(points.size() == 2)
            {
                QLineF l(points.front(), points.back());

                auto midpoint = l.pointAt(0.5);
                auto direction = QVector2D (l.p1() - l.p2()).normalized();
                QTransform rotation; rotation.rotate(90);
                direction = QVector2D(direction.toPointF() * rotation);

                QLineF ray(midpoint, midpoint + (direction.toPointF() * 10000));

                // Convert boundary into segments
                QVector<QLineF> segments;
                auto boundary = polys.front();
                QPointF p1 = boundary.at(0);
                for(int i = 1; i <= boundary.size(); i++){
                    QPointF p2 = boundary.at(i % boundary.size());
                    segments << QLineF(p1,p2);
                    p1 = p2;
                }

                // Intersect with boundary
                QPointF closest;
                double minDist = DBL_MAX;
                for(auto segment : segments)
                {
                    QPointF isect;
                    if(QLineF::BoundedIntersection != ray.intersect(segment, &isect)) continue;
                    double curDist = QVector2D(isect - midpoint).length();
                    if(curDist < minDist){
                        minDist = curDist;
                        closest = isect;
                    }
                }

				// Create smooth curve
				int resolution = 100;
				std::vector<QPointF> points;
				//points.push_back(midpoint + (midpoint - closest));
				points.push_back(l.p1());
				points.push_back(closest);
				points.push_back(l.p2());
				//points.push_back(midpoint + (midpoint - closest));

				SimpleSpline<QPointF> spline(points);
				auto path_points = spline.sampled( resolution );

                QPainterPath loop_path;
				loop_path.moveTo(l.p1());
				for (auto p : path_points) loop_path.lineTo(p);
				//debugs << QVariant::fromValue(loop_path);

				// Draw different results
				auto img = images.front();

				// 1) Loop
				{
					QImage im = img;
					QPainter painter(&im);
					painter.setPen(QPen(Qt::white, 5));
					painter.drawPath(loop_path);

					images << MyQImage(im, 0, img.height());
				}

				// 2) stroke
				{
					QImage im = img;
					QPainter painter(&im);
					painter.setPen(QPen(Qt::white, 5));
					painter.drawLine(l);

					images << MyQImage(im, im.width(), img.height());
				}
				
				// 3) Cut the loop
				{
					QImage im = img;
					QPainter painter(&im);
					painter.setPen(QPen(Qt::white, 5));
					painter.drawPath(loop_path);

					QPainterPath small_boundary;
					small_boundary.addPolygon(polys.back());

					auto closestPointTo = [&](const QPointF &target, const QPainterPath &sourcePath){
						QPointF shortestDistance = sourcePath.elementAt(0) - target;
						qreal shortestLength = shortestDistance.manhattanLength();
						for (int i = 1; i < sourcePath.elementCount(); ++i) {
							const QPointF distance(sourcePath.elementAt(i) - target);
							const qreal length = distance.manhattanLength();
							if (length < shortestLength) {
								shortestDistance = sourcePath.elementAt(i);
								shortestLength = length;
							}
						}
						return shortestDistance;
					};

					auto other_side = closestPointTo(midpoint, small_boundary);
					QLineF loop_line = QLineF(other_side, closest);

					painter.drawLine(loop_line);
					auto delta = loop_line.normalVector();
					delta.setLength(10);
					auto delta_p = delta.p2() - delta.p1();
					
					painter.setPen(QPen(Qt::red, 6));
					painter.drawPoint(loop_line.pointAt(0.5) + delta_p);
					painter.drawPoint(loop_line.pointAt(0.5) - delta_p);

					images << MyQImage(im, im.width() * 2, img.height());
				}

            }
        }

        if(event->modifiers().testFlag(Qt::ShiftModifier))
        {

        }
    }

    update();
}

void Viewer::keyPressEvent(QKeyEvent *event)
{
    if(event->key() == Qt::Key_I) //LOAD IMAGE
    {
        images.clear();
        auto img = QImage(QFileDialog::getOpenFileName(0,"Load image", "", "*.png"));
		auto finalImage = img.width() > 500 ? img.scaledToWidth(500) : img;
		finalImage = finalImage.height() > 500 ? finalImage.scaledToHeight(500) : finalImage;

		// add white pixel border
		QPainter painter(&finalImage);
		painter.setPen(QPen(Qt::white, 3));
		painter.drawRect(finalImage.rect());

		images << MyQImage(finalImage);

        // Collect any possible two boundaries:
        typedef std::vector< std::vector<double> > GridType;
        auto grid = make_grid<double>(images.front());

        typedef std::vector< std::array<double,2> > BoundaryType;
        QVector< BoundaryType > boundaries;
        boundaries << MarchingSquares< std::array<double,2>, GridType >::march(grid, 0);
        floodfill_stacked(0, 0, grid);
        boundaries << MarchingSquares< std::array<double,2>, GridType >::march(grid, 1);

        auto toPoly = [&](const BoundaryType & b){
            QPolygonF poly;
            for(auto coord : b) poly << QPointF(coord[0], coord[1]);
            return poly;
        };

        polys.clear();
        polys << toPoly(boundaries.front());
        if(boundaries.size() > 1) polys << toPoly(boundaries.back());
    }

    if(event->key() == Qt::Key_S)
    {
        if(points.size() < 2) return;

    }

    //Compute Distances toBoundary
    if(event->key() == Qt::Key_D)
    {
    }

    if(event->key() == Qt::Key_U)
    {

    }

    update();
}
