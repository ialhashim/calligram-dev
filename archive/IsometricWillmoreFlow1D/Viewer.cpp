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

    glEnable (GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glEnable( GL_POINT_SMOOTH );
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
    glScaled(0.75,0.75,0.75);

    auto drawPolygon = [&](std::vector<Mesh::Vertex> & vertices, double alpha, QColor c)
    {
        alpha = std::max(0.1, alpha);

        // Lines
        glColor4d(0.4,0.4,0.4,alpha * 0.5);
        glLineWidth(2);
        glBegin(GL_LINE_STRIP);
        for(auto & v : vertices) glVertex3dv(v.data());
        glVertex3dv(vertices.front().data());
        glEnd();

        // Points
        glColor4d(c.redF(),c.greenF(),c.blueF(),alpha);
        glPointSize(4);
        glBegin(GL_POINTS);
        for(auto & v : vertices) glVertex3dv(v.data());
        glEnd();

        // Feature points
        glColor4d(0,0,1,1);
        glPointSize(6);
        glBegin(GL_POINTS);
        glVertex3dv(vertices[0].data());
        glVertex3dv(vertices[vertices.size() * 0.75].data());
        glVertex3dv(vertices[vertices.size() * 0.5].data());
        glVertex3dv(vertices[vertices.size() * 0.25].data());
        glEnd();
    };

    drawPolygon(flow->mesh->flow_vertices[0], 1.0, Qt::yellow);

    for(int i = 0; i < flow->mesh->flow_vertices.size()-1; i++)
    {
        double alpha = double(i) / (flow->mesh->flow_vertices.size()-1);
        drawPolygon(flow->mesh->flow_vertices[i], alpha, QColor(128,0,128));
    }

    drawPolygon(flow->mesh->vertices, 1.0, Qt::red);
}

void Viewer::keyPressEvent(QKeyEvent *e)
{
    if(e->key() == Qt::Key_Space){
        flow->integrate(step_size);
        flow->mesh->center();

        flow->mesh->flow_vertices.push_back(flow->mesh->vertices);

        update();
    }
}
