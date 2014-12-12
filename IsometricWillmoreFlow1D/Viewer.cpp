#include "Viewer.h"

double step_size = 0.01;

Viewer::Viewer(QSharedPointer<IsometricWillmoreFlow1D> flow) : flow(flow)
{
    setMinimumSize(800,800);
}

void Viewer::initializeGL()
{
    initializeOpenGLFunctions();
    glClearColor(1,1,1,1);
}

void Viewer::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Setup 2D
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    // Zoom
    glScaled(0.8,0.8,0.8);

    // Lines
    glColor3d(0.5,0.5,0.5);
    glLineWidth(3);
    glBegin(GL_LINE_STRIP);
    for(auto & v : flow->mesh->vertices) glVertex3dv(v.data());
    glVertex3dv(flow->mesh->vertices.front().data());
    glEnd();

    // Points
    glColor3d(1,0,0);
    glPointSize(6);
    glBegin(GL_POINTS);
    for(auto & v : flow->mesh->vertices) glVertex3dv(v.data());
    glEnd();
}

void Viewer::keyPressEvent(QKeyEvent *e)
{
    if(e->key() == Qt::Key_Space){
        flow->integrate(step_size);
        flow->mesh->center();
        update();
    }
}
