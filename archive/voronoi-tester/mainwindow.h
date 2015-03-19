#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QPainter>
#include <QMap>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private:
    Ui::MainWindow *ui;
};

class PaintWidget : public QWidget {
    void paintEvent(QPaintEvent *);
    void mousePressEvent ( QMouseEvent * event );
    void keyPressEvent( QKeyEvent * event );
    void wheelEvent(QWheelEvent *event);
public:
    PaintWidget(QWidget * parent = 0);
    QPolygonF points, innerPoints;
    QString defaultFile;
    QVector<QPainterPath> paths;
    QMap<QString,QPolygonF> polygons;
    QVector<QLineF> lineSegments;
    QVector<QLineF> closest;

	QVector<QLineF> triangulation;

    void computeVoronoi();
};

#endif // MAINWINDOW_H
