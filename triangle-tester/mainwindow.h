#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

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
public:
    PaintWidget(QWidget * parent = 0);

	QVector<QPolygonF> triangles;
};

#endif // MAINWINDOW_H
