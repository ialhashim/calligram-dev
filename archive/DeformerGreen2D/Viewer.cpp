#include "Viewer.h"
#include "ui_Viewer.h"
#include <QPainter>
#include <QDebug>
#include <QTimer>
#include <QKeyEvent>
#include <QFileDialog>

#include "mvc.h"

#include "Mesh.h"
using namespace surface_mesh;

inline Eigen::Vector3d toEigen(QPointF & p){ return Eigen::Vector3d(p.x(), p.y(), 0); }

Viewer::Viewer(QWidget *parent) : QGLWidget(parent), ui(new Ui::Viewer), isCoordinatesReady(false)
{
	QGLFormat glf = QGLFormat::defaultFormat();
	glf.setSampleBuffers(true);
	this->setFormat(glf);

    ui->setupUi(this);
	setFocusPolicy(Qt::ClickFocus);

    // Experiment
    demo();
}

void Viewer::demo()
{
    meshes.push_back(new Mesh("C:/Development/calligram-dev/letterMeshes/A.obj", 0.3, Point(100,100,0)));
}

Viewer::~Viewer()
{
    delete ui;
}

void Viewer::paintEvent(QPaintEvent *)
{
    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing);
    painter.fillRect(rect(),Qt::white);

	// CHULL
	QPainterPath path2;
	path2.addPolygon(chull);
	painter.setPen(QPen(Qt::blue, 2));
	painter.drawPath(path2);
	painter.setPen(QPen(Qt::green, 8));
	painter.drawPoints(chull);

	// Draw points
	painter.setPen(QPen(Qt::black, 1));
    QPainterPath path;
    path.addPolygon(points);
    painter.drawPath( path );
    painter.setPen(QPen(Qt::red, 6));
    painter.drawPoints( points );

    // Draw mesh
    {
        painter.beginNativePainting();

        for(Mesh * m : meshes)
        {
            auto & mesh = *m;

            //glViewport(pos.x(), pos.y() - height, width, height);
            glMatrixMode(GL_PROJECTION);
            glPushMatrix();
            glLoadIdentity();
            glOrtho(0, this->width(), this->height(), 0, 0.0, -1.0);
            glMatrixMode(GL_MODELVIEW);
            glPushMatrix();
            glLoadIdentity();

            glColor3d(1,0,0);

            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            glBegin(GL_TRIANGLES);

            auto & points = mesh.vertex_property<Point>("v:point");
            for(auto triangle : mesh.faces())
            {
                std::vector<Point> tripnts;

                for(auto v : Surface_mesh::Vertex_around_face_circulator(m, triangle))
                    tripnts.push_back(points[v]);

                glVertex3dv( tripnts[0].data() );
                glVertex3dv( tripnts[1].data() );
                glVertex3dv( tripnts[2].data() );
            }
            glEnd();
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        }

        painter.endNativePainting();
    }
}

void Viewer::mousePressEvent(QMouseEvent *event)
{
    if(event->buttons().testFlag(Qt::LeftButton))
    {
        if(event->modifiers().testFlag(Qt::ControlModifier)){
            if(!points.empty()) points.pop_back();
            points << event->pos();
            points << points.front();
        }

        if(event->modifiers().testFlag(Qt::ShiftModifier)){
            points.clear();
        }
    }

    update();
}

void Viewer::keyPressEvent(QKeyEvent *event)
{
	if(event->key() == Qt::Key_Space){
		QTimer * timer = new QTimer;
		connect(timer, &QTimer::timeout, [=]()
		{
			if (!isCoordinatesReady) computeGreenCoordinates();

			static QPolygonF origPositions = points;
			static double theta = 0.0;
			static double radius = 50.0;

			// Move some control point
			for (int i = 0; i < 3; i++)	points[i] = origPositions[i] + (QPointF(cos(theta), sin(theta)) * radius);
			
			applyDeformGreen();

			update();

			theta += 0.1;
		});
		timer->start( 10 );
	}

	if (event->key() == Qt::Key_Backspace){
		QTimer * timer = new QTimer;

		// Compute MVC weights:
		if (meshes.empty() || points.empty()) return;
		Mesh & mesh = *meshes.front();
		auto & mesh_points = mesh.vertex_property<Point>("v:point");
		auto & mesh_weights = mesh.vertex_property<std::vector<double> >("v:mvc", std::vector<double>());

		std::vector<Eigen::Vector2d> cage;
		for (auto p : points) cage.push_back(Eigen::Vector2d(p.x(), p.y()));

		for (auto v : mesh.vertices()){
			auto p = mesh_points[v];
			mesh_weights[v] = MeanValueCoordinates<Eigen::Vector2d>::computeWeights(p[0], p[1], cage);
		}

		connect(timer, &QTimer::timeout, [=]()
		{
			static QPolygonF origPositions = points;
			static double theta = 0.0;
			static double radius = 50.0;

			// Move some control point
			for (int i = 0; i < 3; i++)	points[i] = origPositions[i] + (QPointF(cos(theta), sin(theta)) * radius);

			applyDeformMVC();

			update();
			theta += 0.1;
		});
		timer->start(10);
	}

	// Compute GC coordinates for current cage
	if (event->key() == Qt::Key_C)
	{
		computeGreenCoordinates();
	}

	QString DELIM = " ";

    if(event->key() == Qt::Key_O){
        // Load points from file
        QFile file( QFileDialog::getOpenFileName(0, "Load points", "","Points File (*.pts)") );
        file.open(QFile::ReadOnly | QFile::Text);
        QTextStream in(&file);
        QStringList lines = in.readAll().split('\n');
        for( auto line : lines ){
            auto pointCoord = line.split(DELIM);
            if(pointCoord.size() < 2) continue;
            points << QPointF( pointCoord[0].toDouble(), pointCoord[1].toDouble() );
        }
    }

    if(event->key() == Qt::Key_S){
        if(points.empty()) return;

        QFile file( QFileDialog::getSaveFileName(0, "Save points", "","Points File (*.pts)") );
        file.open(QFile::WriteOnly | QFile::Text);
        QTextStream out(&file);
        for(auto p : points){
			out << p.x() << DELIM << p.y() << "\n";
        }
    }

    // Test convexhull:
    if(event->key() == Qt::Key_H){
        std::vector<Eigen::Vector3d> convexhull;

        std::vector<Eigen::Vector3d> pnts;
        for(auto p : points) pnts.push_back(toEigen(p));

		convexhull2d<Eigen::Vector3d>(pnts.begin(), pnts.end(), std::back_inserter(convexhull));

		for (auto p : convexhull) chull << QPointF(p[0],p[1]);
    }

    update();
}

void Viewer::computeGreenCoordinates()
{
	if (meshes.empty() || points.empty()) return;
	Mesh & mesh = *meshes.front();

	auto & mesh_points = mesh.vertex_property<Point>("v:point");
	auto & gcoord_phi = mesh.vertex_property<std::vector<double> >("v:green_coord_phi", std::vector<double>());
	auto & gcoord_psi = mesh.vertex_property<std::vector<double> >("v:green_coord_psi", std::vector<double>());

	int n_boundary_pt = points.size();

	/// Compute properties for cage:
	{
		// Orientation:
		double area = 0;
		for (int i = 0; i < n_boundary_pt; ++i)
		{
			const Vector3 &v1 = toEigen(points[i]); // current
			const Vector3 &v2 = toEigen(points[(i + 1) % n_boundary_pt]); // next
			const Vector3 &v3 = toEigen(points[(i - 1 + n_boundary_pt) % n_boundary_pt]); // previous
			area += v1[0] * (v2[1] - v3[1]);
			area *= 0.5;
		}
		if (area < 0)
		{
			std::reverse(points.begin(), points.end());
		}

		// Cage edge properties:
		control_edge_length.resize(n_boundary_pt);
		control_edge_normal.resize(n_boundary_pt);

		for (int i = 0; i < n_boundary_pt; ++i)
		{
			const Vector3 &p0 = toEigen(points[i]); // current
			const Vector3 &p1 = toEigen(points[(i + 1) % n_boundary_pt]); // next
			control_edge_length[i] = (p0 - p1).norm();

			double x = p1[0] - p0[0];
			double y = p1[1] - p0[1];

			double len = sqrt(x*x + y*y);
			control_edge_normal[i] = QPointF(y / len, -x / len);
		}
	}

	/// Compute weights for each mesh vertex:
	for (auto vertex : mesh.vertices())
	{
		Vector3 queryCoord = mesh_points[vertex];

		Eigen::Vector2d a, b, nj;
		double V, Q, S, R, BA, SRT, L0, L1, A0, A1, A10, L10;

		std::vector<double> phi(n_boundary_pt, 0.0);
		std::vector<double> psi(n_boundary_pt, 0.0);

		for (int j = 0; j < n_boundary_pt; j++)
		{
			int jp = (j + 1) % n_boundary_pt;
			const Vector3 &vj1 = toEigen(points[j]);
			const Vector3 &vj2 = toEigen(points[jp]);

			a[0] = vj2[0] - vj1[0];
			a[1] = vj2[1] - vj1[1];

			nj[0] = a[1];
			nj[1] = -a[0];

			b[0] = vj1[0] - queryCoord[0];
			b[1] = vj1[1] - queryCoord[1];

			Q = a[0] * a[0] + a[1] * a[1];
			S = b[0] * b[0] + b[1] * b[1];
			R = 2.0*(a[0] * b[0] + a[1] * b[1]);
			BA = b[0] * nj[0] + b[1] * nj[1];   // Remember do not normalize Normal vector here.
			V = 4.0*S*Q - R*R;
			assert(V > 0.0);
			SRT = sqrt(V);
			assert(SRT > 0.0);
			L0 = log(S);
			L1 = log(S + Q + R);
			A0 = atan(R / SRT) / SRT;
			A1 = atan((2 * Q + R) / SRT) / SRT;
			A10 = A1 - A0;
			L10 = L1 - L0;
			psi[j] = sqrt(Q) / (4.0*M_PI)*((4.0*S - (R*R / Q))*A10 + (R / (2.0*Q))*L10 + L1 - 2);
			phi[jp] -= (BA / (2.0*M_PI))*((L10 / (2.0*Q)) - A10*R / Q);
			phi[j] += (BA / (2.0*M_PI))*((L10 / (2.0*Q)) - A10*(2.0 + R / Q));
		}

		double sum = 0.0;
		for (int i = 0; i < n_boundary_pt; ++i)
			sum += phi[i];

		for (int i = 0; i < n_boundary_pt; ++i)
			phi[i] /= sum;

		sum = 0.0;
		for (int i = 0; i < n_boundary_pt; ++i)
			sum += phi[i];
		assert(fabs(sum - 1.0) < 1.0E-06);

		gcoord_phi[vertex] = phi;
		gcoord_psi[vertex] = psi;
	}

	isCoordinatesReady = true;
}

void Viewer::applyDeformGreen()
{
	if (meshes.empty() || points.empty()) return;
	if (!isCoordinatesReady) computeGreenCoordinates();
	
	Mesh & mesh = *meshes.front();
	auto & mesh_points = mesh.vertex_property<Point>("v:point");
	auto & gcoord_phi = mesh.vertex_property< std::vector<double> >("v:green_coord_phi");
	auto & gcoord_psi = mesh.vertex_property< std::vector<double> >("v:green_coord_psi");

	int n_boundary_pt = points.size();

	auto cur_control_edge_length = control_edge_length;
	for (int i = 0; i < n_boundary_pt; ++i)
	{
		const Vector3 &p0 = toEigen(points[i]); // current
		const Vector3 &p1 = toEigen(points[(i + 1) % n_boundary_pt]); // next
		cur_control_edge_length[i] = (p0 - p1).norm();
	}

	for (auto vertex : mesh.vertices())
	{
		Vector3 weighted_pt(0, 0, 0);

		for (int j = 0; j < n_boundary_pt; j++)
		{
			auto cagepnt_pos = toEigen(points[j]);
			auto phi = gcoord_phi[vertex][j]; 
			auto psi = gcoord_psi[vertex][j];

			//OMesh::Point pos = point_data[controlpt_idx_to_vtx_idx_map[j]];
			//weighted_pt += greenw(i, j)*pos + greenw(i, j + n_boundary_pt)*cur_control_edge_length[j] * control_edge_normal[j] / control_edge_length[j];

			weighted_pt += (phi * cagepnt_pos) + (psi * cur_control_edge_length[j] * toEigen(control_edge_normal[j]) / control_edge_length[j]);
		}

		mesh_points[vertex] = weighted_pt;
	}
}

void Viewer::applyDeformMVC()
{
	if (meshes.empty() || points.empty()) return;
	Mesh & mesh = *meshes.front();
	auto & mesh_points = mesh.vertex_property<Point>("v:point");
	auto & mesh_weights = mesh.vertex_property<std::vector<double> >("v:mvc");

	std::vector<Eigen::Vector2d> cage;
	for (auto p : points) cage.push_back(Eigen::Vector2d(p.x(), p.y()));

	for (auto v : mesh.vertices())
	{
		auto pp = MeanValueCoordinates<Eigen::Vector2d>::interpolate2d(mesh_weights[v], cage);
		mesh_points[v] = Eigen::Vector3d(pp[0],pp[1],0);
	}
}
