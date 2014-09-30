#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QMouseEvent>
#include <QFileInfo>
#include <QTextStream>

#include "voronoi.h"

MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->mainlayout->addWidget( new PaintWidget(this) );
}

MainWindow::~MainWindow()
{
    delete ui;
}

void PaintWidget::paintEvent(QPaintEvent *) {
    QPainter p(this);
    p.setRenderHint(QPainter::Antialiasing);
    p.setRenderHint(QPainter::HighQualityAntialiasing);

    QPainterPath path;
    path.addPolygon(points);
    p.drawPath( path );
    p.setPen(QPen(Qt::red, 3));
    p.drawPoints( points );

    // Voronoi
    p.setPen(QPen(Qt::blue, 5));
    p.drawPoints( polygons["voroVerts"] );

    // Line segments
    p.setPen(QPen(Qt::black, 1));
    for(auto l : lineSegments){
        p.drawLine(l);
    }
}

void PaintWidget::mousePressEvent(QMouseEvent *event)
{
    if(event->buttons().testFlag(Qt::LeftButton))
    {
        if(event->modifiers().testFlag(Qt::ControlModifier)){
            if(!points.empty()) points.pop_back();
            points << event->pos();
            points << points.front();
        }

        if(event->modifiers().testFlag(Qt::ShiftModifier)){
            points.clear();
        }
    }
    else
    {
        QFile file( defaultFile );
        file.open(QFile::WriteOnly | QFile::Text);
        QTextStream out(&file);
        QPolygonF copyPoints = points;
        copyPoints.removeLast();
        for(auto p : copyPoints){
            out << p.x() << "," << p.y() << "\n";
        }
    }

    update();
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

void PaintWidget::keyPressEvent(QKeyEvent *event)
{
    if(event->key() == Qt::Key_Space)
    {
        // Resample path
        points = resamplePolygon(points, 300);
    }

    if(event->key() == Qt::Key_S)
    {
        // Smooth path
        points = smoothPolygon(points, 4);
    }

    if(event->key() == Qt::Key_V)
    {
        computeVoronoi();
    }

    update();
}

void PaintWidget::wheelEvent(QWheelEvent *event)
{
    static double scaleFactor = 1.0;

    int numSteps = event->delta() / 15 / 8;
    scaleFactor = pow(1.25, numSteps);

    QTransform t;
    t.scale( scaleFactor, scaleFactor );
    points = t.map(points);

    update();
}

PaintWidget::PaintWidget(QWidget *parent) : QWidget(parent) {
    defaultFile = "points.txt";

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

    // Options
    setFocusPolicy(Qt::ClickFocus);
}

void PaintWidget::computeVoronoi()
{
    QPolygonF pointsCopy = points;
    pointsCopy.removeLast();

    // Format values
    QVector< QVector<float> > vals(2);
    for(auto p : pointsCopy){
        vals[0].push_back( p.x() );
        vals[1].push_back( p.y() );
    }
    QSizeF pointsSize = points.boundingRect().size();
    float bound = std::max( pointsSize.width(), pointsSize.height() ) * 2;

    // Voronoi diagram
    VoronoiDiagramGenerator vdg;
    vdg.generateVoronoi(&vals[0][0], &vals[1][0], vals[0].size(), -bound, bound, -bound, bound, 1e-6f);

    QPolygonF voroVerts;
    float x,y;
    while(vdg.getNextVertex(x,y)){
        voroVerts.push_back( QPointF(x,y) );
    }
    polygons["voroVerts"] = voroVerts;

    lineSegments.clear();
    float x1,y1,x2,y2;
    vdg.resetIterator();
    while(vdg.getNext(x1,y1,x2,y2)){
        QPointF p1(x1,y1), p2(x2,y2);

        // Only inner edges
        if( points.containsPoint(p1, Qt::WindingFill) && points.containsPoint(p2, Qt::WindingFill) ){
            lineSegments.push_back( QLineF(p1,p2) );
        }
    }
}
