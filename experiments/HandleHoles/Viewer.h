#pragma once

#include <QGLWidget>
#include <QVector>

namespace Ui {
class Viewer;
}

class MyQImage
        : public QImage{public: int x,y; double opacity; MyQImage(const QImage & img = QImage(),
                               int x = 0, int y = 0, double opacity = 1.0) : QImage(img), x(x), y(y), opacity(opacity) {}};

class Viewer : public QGLWidget
{
    Q_OBJECT

public:
    explicit Viewer(QWidget *parent = 0);
    ~Viewer();
    void paintEvent(QPaintEvent *);
    void keyPressEvent(QKeyEvent *);
    void mousePressEvent(QMouseEvent *);

    QVector<MyQImage> images;
    QPolygonF points;

    QVector<QPolygonF> polys;

    QVector<QVariant> debugs;

private:
    Ui::Viewer *ui;

public slots:
};
