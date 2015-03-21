#include "Viewer.h"
#include "Mesh.h"

#include <QTimer>

#include "globals.h"
#include "trianglelib.h"
#include "libShapeOp/Solver.h"
#include "libShapeOp/Force.h"
#include "libShapeOp/Constraint.h"

QSharedPointer<surface_mesh::Mesh> _mesh;

std::vector<int> pins_const_ids;
std::vector<int> pins;

Viewer::Viewer()
{
    setMinimumSize(800,800);
    setFocusPolicy(Qt::ClickFocus);
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

	glScaled(1, 1, 1); // Zoom
	glTranslated(-0.5, -0.2, 0);

	glPointSize(5);
	glLineWidth(5);
	glColor3d(0, 0, 0);

    // Draw mesh
    //glPolygonMode(GL_FRONT, GL_LINE);
    if(!mesh.isNull() && mesh->n_faces())
    {
        Vector3 sum(0,0,0);
        for(auto pi : mesh->vertices())
        {
            auto blah = mesh->get_vertex_property<Vector3>("v:point")[pi];
            sum += blah;
        }
        Vector3 c = sum / mesh->n_vertices();

        double delta = -0.5;
        c += Vector3(delta, delta,0);

        auto points = mesh->vertex_property<Eigen::Vector3d>("v:point");

        glBegin(GL_TRIANGLES);
        for(auto f : mesh->faces())
        {
            for(auto v : surface_mesh::Mesh::Vertex_around_face_circulator(mesh.data(), f))
            {
                auto realPos = points[v];

                auto fixedPos = realPos - c;

                glVertex2d(fixedPos[0], fixedPos[1]);
            }
        }
        glEnd();

        glColor3d(1,0,0);
        glBegin(GL_POINTS);
        for(auto idx : pins)
        {
            auto p = points[surface_mesh::Surface_mesh::Vertex(idx)];
            p -= c;
            glVertex2d(p[0], p[1]);
        }
        glEnd();
    }
}

void Viewer::keyPressEvent(QKeyEvent *e)
{
	if (e->key() == Qt::Key_S){
		this->grabFramebuffer().save("captured.png");
	}

    if(e->key() == Qt::Key_Space){
        trianglelib::BasicMesh<double> m;
        LetterGeometry letter('p', 70, 0.02);
        trianglelib::triangulatePolygon<double>(letter.points, letter.holes, m);

        mesh = QSharedPointer<surface_mesh::Mesh>( new surface_mesh::Mesh() );

        for(auto v : m.vertices) mesh->add_vertex(Eigen::Vector3d(v[0],v[1],0));
		for (auto t : m.triangles) mesh->add_triangle(t[0], t[1], t[2]);

		// Solver
		{
			int dimensions = 2;

			solver = QSharedPointer<ShapeOp::Solver>(new ShapeOp::Solver());
			size_t nb_points = mesh->n_vertices();
			auto points = mesh->vertex_property<Eigen::Vector3d>("v:point");
			Eigen::Map<ShapeOp::Matrix3X> p(points.vector().front().data(), 3, nb_points);

			solver->setPoints(p);

			{
				double triangle_weight = 0.2;

				for (auto & face : mesh->faces())
				{
					std::vector<int> id_vector;
					for (auto & v : mesh->vertices(face)) id_vector.push_back(v.idx());

					auto c = std::make_shared<ShapeOp::TriangleStrainConstraint>(id_vector, triangle_weight, p, dimensions == 2);
					solver->addConstraint(c);
				}
			}
			{
				double bending_weight = 0;

				for (auto & edge : mesh->edges())
				{
					if (mesh->is_boundary(edge)) continue;

					auto h2 = mesh->halfedge(edge, 1);
					auto h3 = mesh->halfedge(edge, 0);

					auto v2 = mesh->to_vertex(h2);
					auto v3 = mesh->to_vertex(h3);
					auto v1 = mesh->to_vertex(mesh->next_halfedge(h3));
					auto v4 = mesh->to_vertex(mesh->next_halfedge(h2));

					std::vector<int> id_vector;
					id_vector.push_back(v2.idx());
					id_vector.push_back(v3.idx());
					id_vector.push_back(v1.idx());
					id_vector.push_back(v4.idx());

					auto c = std::make_shared<ShapeOp::BendingConstraint>(id_vector, bending_weight, p);
					solver->addConstraint(c);
				}
			}
			{
				double area_weight = 0;

				for (auto & face : mesh->faces())
				{
					std::vector<int> id_vector;
					for (auto & v : mesh->vertices(face)) id_vector.push_back(v.idx());

					auto c = std::make_shared<ShapeOp::AreaConstraint>(id_vector, area_weight, p);
					solver->addConstraint(c);
				}
			}
			{
				double laplacian_weight = 0;

				for (auto & v : mesh->vertices())
				{
					//if (mesh->is_boundary(v)) continue;

					std::vector<int> id_vector;
					id_vector.push_back(v.idx());
					for (auto & h : surface_mesh::Mesh::Halfedge_around_vertex_circulator(mesh.data(), v)) 
						id_vector.push_back(mesh->to_vertex(h).idx());

					auto c = std::make_shared<ShapeOp::UniformLaplacianConstraint>(id_vector, laplacian_weight, p, false);
					solver->addConstraint(c);
				}
			}

			// Constraints

			pins.push_back(0);
            pins.push_back(mesh->n_vertices() * 0.5);
            pins.push_back(mesh->n_vertices() * 0.25);
            pins.push_back(mesh->n_vertices() * 0.75);
            pins.push_back(mesh->n_vertices() * 0.4);

			{
				double close_weight = 1.0;

				for (auto pinned : pins)
				{
					std::vector<int> id_vector;
					id_vector.push_back(pinned);
					auto c = std::make_shared<ShapeOp::ClosenessConstraint>(id_vector, close_weight, p);
					auto cid = solver->addConstraint(c);

                    pins_const_ids.push_back(cid);
				}
			}

			// Gravity
			auto f = std::make_shared<ShapeOp::GravityForce>(Vector3(0, -0.01, 0));
            //solver->addForces(f);

			/// 4) Initialize the solver
			solver->initialize(true, 0.01, 0.5, 1.0);

			/// 5) Optimize
			solver->solve(1, dimensions);

			// Solution
			auto final_points = solver->getPoints();
			for (size_t i = 0; i < p.cols(); i++) p.col(i) = final_points.col(i);

			auto timer = new QTimer(this);
			connect(timer, &QTimer::timeout, [&]() {
				size_t nb_points = this->mesh->n_vertices();
				auto points = mesh->vertex_property<Eigen::Vector3d>("v:point");
				Eigen::Map<ShapeOp::Matrix3X> p(points.vector().front().data(), 3, nb_points);
				int dimensions = 2;

                solver->solve(20, dimensions);


				auto final_points = solver->getPoints();
				for (size_t i = 0; i < p.cols(); i++) p.col(i) = final_points.col(i);
				update();

                for(int i = 0; i < pins.size(); i++)
                {
                    // Modifty constraints
                    int idx = i;
                    auto c = std::dynamic_pointer_cast <ShapeOp::ClosenessConstraint>(solver->getConstraint(pins_const_ids[idx]));

                    auto pts = solver->getPoints();
                    Vector3 oldPos = pts.col(pins[idx]);
                    Vector3 newPos = oldPos + (Vector3::Random() * 0.2 * (double(rand())/RAND_MAX));
                    c->setPosition( newPos );

                    //if((double(rand())/RAND_MAX) > 0.5) continue;
                }

			});
			timer->start(30);
		}

        update();
    }
}
