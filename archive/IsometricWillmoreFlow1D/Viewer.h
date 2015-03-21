#ifndef VIEWER_H
#define VIEWER_H

#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QSharedPointer>
#include <QKeyEvent>
#include "IsometricWillmoreFlow1D.h"

class Viewer : public QOpenGLWidget, protected QOpenGLFunctions
{
    Q_OBJECT
public:
    Viewer(QSharedPointer<IsometricWillmoreFlow1D> flow);
    QSharedPointer<IsometricWillmoreFlow1D> flow;

protected:
    void initializeGL();
    void paintGL();

    void keyPressEvent(QKeyEvent*e);
};

#endif // VIEWER_H
