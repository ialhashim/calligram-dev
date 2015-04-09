#include "Viewer.h"
#include "ui_Viewer.h"
#include <QPainter>
#include <QDebug>
#include <QTimer>
#include <QKeyEvent>
#include <QMouseEvent>
#include <QFileDialog>

#include "libfastmarching.h"
libfastmarching::DynamicGrid * dgrid = nullptr;

#include "globals.h"

Viewer::Viewer(QWidget *parent) : QGLWidget(parent), ui(new Ui::Viewer)
{
	QGLFormat glf = QGLFormat::defaultFormat();
	glf.setSampleBuffers(true);
	this->setFormat(glf);

    ui->setupUi(this);
	setFocusPolicy(Qt::ClickFocus);

    QString defaultImage = "n_outline.png";
    if(QFileInfo(defaultImage).exists()) loadImage(defaultImage);
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

void Viewer::mousePressEvent(QMouseEvent *event)
{
    if(event->buttons().testFlag(Qt::LeftButton))
    {
        if(event->modifiers().testFlag(Qt::ControlModifier))
        {
            auto shape = images.front();

            if(qRed(shape.pixel(event->pos()))) points << event->pos();
        }

        if(event->modifiers().testFlag(Qt::ShiftModifier))
        {
            points.clear();
        }
    }

    update();
}

void Viewer::loadImage(QString filename)
{
    auto img = QImage(filename);

    images.clear();
    images << MyQImage(img);
    points.clear();

    if(dgrid) delete dgrid;
    dgrid = new libfastmarching::DynamicGrid(img);
}

void Viewer::keyPressEvent(QKeyEvent *event)
{
    if(event->key() == Qt::Key_I)
    {
        this->loadImage(QFileDialog::getOpenFileName(0,"Load image", "", "*.png"));
    }

    if(event->key() == Qt::Key_W)
    {
        auto start_point = points.at(points.size()-2).toPoint();
        auto end_point = points.back().toPoint();

        QImage viz;
        libfastmarching::Path fm2path = dgrid->bestPath(start_point, end_point, viz);

        dgrid->modifyWalls(fm2path, globalStrokeSize * 0.5);

        // Post-process path
        QPolygon path;
        for(auto p : fm2path) path << QPoint(p[0], p[1]);
        path = smoothPolygon(resamplePolygon(path, dgrid->walls.width() * 0.25), 2);

        // Visualize
        images << visualizeStroke(dgrid->walls, path);
    }

    if(event->key() == Qt::Key_S)
    {
        if(points.size() < 2) return;

        auto shape = images.front();
        auto start_point = points.at(points.size()-2).toPoint();
        auto end_point = points.back().toPoint();

        libfastmarching::Path fm2path;
        libfastmarching::fm2star(shape, start_point, end_point, fm2path);

        // Post-process path
        QPolygon path;
        for(auto p : fm2path) path << QPoint(p[0], p[1]);
        path = smoothPolygon(resamplePolygon(path,shape.width() * 0.25), 2);

        // Visualize
        images << visualizeStroke(shape, path);
    }

    update();
}
