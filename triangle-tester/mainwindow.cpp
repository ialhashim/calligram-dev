#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QPainter>
#include <QFileInfo>
#include <QTextStream>
#include "iostream"

#include "trianglelib.h"

using namespace std;

MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->mainlayout->addWidget( new PaintWidget(this) );
}

MainWindow::~MainWindow()
{
    delete ui;
}

PaintWidget::PaintWidget(QWidget * parent) : QWidget(parent)
{
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

    std::vector< std::vector<double> > all_points;

    for(auto p : points){
        std::vector<double> point(2, 0);
        point[0] = p.x();
        point[1] = p.y();
        all_points.push_back(point);
    }

    trianglelib::BasicMesh<double> mesh;
    trianglelib::triangulatePolygon<double>(all_points, mesh);

	for(auto triangle : mesh.triangles)
	{
		QPolygonF tri;

		tri << QPointF( mesh.vertices[triangle[0]][0], mesh.vertices[triangle[0]][1] );
		tri << QPointF( mesh.vertices[triangle[1]][0], mesh.vertices[triangle[1]][1] );
		tri << QPointF( mesh.vertices[triangle[2]][0], mesh.vertices[triangle[2]][1] );

        triangles.push_back( tri );
	}
}

void PaintWidget::paintEvent(QPaintEvent *)
{
    QPainter p(this);

	for(auto triangle : triangles)
		p.drawPolygon(triangle);
}
