#pragma once

#include <QGLWidget>
#include <QVector>

namespace surface_mesh{class Mesh;}

namespace Ui {
class Viewer;
}

class MyQImage : public QImage{public: int x,y; double opacity; MyQImage(const QImage & img = QImage(),
                               int x = 0, int y = 0, double opacity = 1.0) : QImage(img), x(x), y(y), opacity(opacity) {}};

class Viewer : public QGLWidget
{
    Q_OBJECT

public:
    explicit Viewer(QWidget *parent = 0);
    ~Viewer();

    void paintEvent(QPaintEvent *);
    void keyPressEvent(QKeyEvent *event);
    void mousePressEvent(QMouseEvent *event);

    QVector<MyQImage> images;

	QVector<surface_mesh::Mesh*> meshes;	   
	
	// Cage
	QPolygonF points;
	QVector<double> control_edge_length;
	QVector<QPointF> control_edge_normal;
	bool isCoordinatesReady;

private:
    Ui::Viewer *ui;

public slots:
    void demo();

	// code from: https://github.com/Juyong/OpenFlipper/Plugin-LBC
	void computeGreenCoordinates();
	void applyDeform();
};
