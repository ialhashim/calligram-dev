#pragma once
#include <QPainter>
#include <QPolygon>
#include <QLine>

#include "trianglelib.h"

struct ARAPDeformer;

class MyPolygon : public QPolygonF
{
public:
    MyPolygon( QVector<QPointF> fromPoints = QVector<QPointF>(), QString id = "polygon" );

    void computeSkeleton();
    void computeMesh();

	void prepareDeform();
	void deform();
	ARAPDeformer * deformer;

    std::vector<QLineF> skeleton;
    trianglelib::BasicMesh<double> mesh;

	void draw(QPainter &painter, QPoint pos, int width, int height);

	QString id;
};
