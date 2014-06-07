#include "Polygon.h"
#include "voronoi.h"
#include <QGLWidget>

#include "ARAPDeformer.h"

MyPolygon::MyPolygon(QVector<QPointF> fromPoints, QString id) : QPolygonF( fromPoints ), id(id), deformer(NULL)
{
    if( fromPoints.size() > 3 )
    {
        computeSkeleton();
        computeMesh();

		prepareDeform();
    }
}

void MyPolygon::computeSkeleton()
{
    QPolygonF points;
    for(auto p : *this) points << p;
    points.removeLast();

    // Format values
    QVector< QVector<float> > vals(2);
    for(auto p : points){
        vals[0].push_back( p.x() );
        vals[1].push_back( p.y() );
    }
    QSizeF pointsSize = points.boundingRect().size();
    float bound = std::max( pointsSize.width(), pointsSize.height() ) * 2;

    // Voronoi diagram
    VoronoiDiagramGenerator vdg;
    vdg.generateVoronoi(&vals[0][0], &vals[1][0], vals[0].size(), -bound, bound, -bound, bound, 1e-6f);

    // Vertices
    {
        QPolygonF voroVerts;
        float x,y;
        while(vdg.getNextVertex(x,y)){
            voroVerts.push_back( QPointF(x,y) );
        }
    }

    // Edges
    {
        skeleton.clear();
        float x1,y1,x2,y2;
        vdg.resetIterator();
        while(vdg.getNext(x1,y1,x2,y2)){
            QPointF p1(x1,y1), p2(x2,y2);

            // Only inner edges
            if( points.containsPoint(p1, Qt::WindingFill) && points.containsPoint(p2, Qt::WindingFill) ){
                skeleton.push_back( QLineF(p1,p2) );
            }
        }
    }
}

void MyPolygon::computeMesh()
{
	QPolygonF points;
	for(auto p : *this) points << p;
	points.removeLast();

    std::vector< std::vector<double> > all_points;
    for(auto p : points){
        std::vector<double> point(2, 0);
        point[0] = p.x();
        point[1] = p.y();
        all_points.push_back(point);
    }

    trianglelib::triangulatePolygon<double>(all_points, mesh, "pq30");
}

void MyPolygon::draw(QPainter &painter, QPoint pos, int width, int height)
{
    // Draw skeleton
	{
		painter.setPen( QPen(Qt::blue,2) );
        for(auto l : skeleton) painter.drawLine(l);
    }

	// Draw mesh
	{
		painter.beginNativePainting();

		glViewport(pos.x(), pos.y() - height, width, height);
		glMatrixMode(GL_PROJECTION);
		glPushMatrix();
		glLoadIdentity();
		glOrtho(0, width, height, 0, 0.0, -1.0);
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glLoadIdentity();

		glColor3d(1,0,0);

		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glBegin(GL_TRIANGLES);
		for(auto triangle : mesh.triangles){
			std::vector<double> v0 = mesh.vertices[triangle[0]];
			std::vector<double> v1 = mesh.vertices[triangle[1]];
			std::vector<double> v2 = mesh.vertices[triangle[2]];

			glVertex2dv( &v0[0] );
			glVertex2dv( &v1[0] );
			glVertex2dv( &v2[0] );
		}
		glEnd();
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

		painter.endNativePainting();
	}
}

void MyPolygon::prepareDeform()
{
	Surface_mesh * m = new Surface_mesh;

	// Add vertices
	for(auto v : mesh.vertices)
	{
		Vector3d point(v[0], v[1], 0.0);
		m->add_vertex( point );
	}

	// Add faces
	for(auto face : mesh.triangles)
	{
		m->add_triangle(Surface_mesh::Vertex(face[0]), Surface_mesh::Vertex(face[1]), Surface_mesh::Vertex(face[2]));
	}

	m->update_face_normals();
	m->update_vertex_normals();

	deformer = new ARAPDeformer(m);
	
	// Set fixed and movable vertices
	deformer->SetAnchor( Surface_mesh::Vertex(m->n_vertices() * 0.5) );
	deformer->SetControl( Surface_mesh::Vertex(0) );
	deformer->Deform(0);

	//deformer->PushControl( Surface_mesh::Vertex(0), Vector3d(-10,-10,0) );
	//deformer->Deform(2);

	// DEBUG:
	m->write( id.toStdString() + ".obj" );
}

void MyPolygon::deform()
{
	static Vector3d origPos = Vector3d(mesh.vertices[0][0], mesh.vertices[0][1], 0.0);
	static double theta = 0.0;
	static double radius = 20.0;

	// Move control point
	deformer->SetControlPosition(Surface_mesh::Vertex(0), origPos + (Vector3d(cos(theta), sin(theta),0) * 10));
	deformer->Deform( 1 );

	theta += 0.1;

	// Update basic mesh
	for(auto v : deformer->mesh->vertices())
	{
		Vector3d newPos = deformer->points[v];
		mesh.vertices[v.idx()][0] = newPos[0];
		mesh.vertices[v.idx()][1] = newPos[1];
	}
}
