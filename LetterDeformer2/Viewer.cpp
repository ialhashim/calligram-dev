#include "Viewer.h"
#include "Mesh.h"

#include <QTimer>
#include <QElapsedTimer>

#include "globals.h"
#include "trianglelib.h"
#include "libShapeOp/Solver.h"
#include "libShapeOp/Force.h"
#include "libShapeOp/Constraint.h"

#include "GeoHeatHelper.h"

QSharedPointer<surface_mesh::Mesh> _mesh;

std::vector<int> pins_const_ids;
std::vector<int> pins;

QTimer * timer = nullptr;
QElapsedTimer time;

Viewer::Viewer()
{
	setMinimumSize(800, 800);
	setFocusPolicy(Qt::ClickFocus);
}

void Viewer::initializeGL()
{
	initializeOpenGLFunctions();
	glClearColor(1, 1, 1, 1);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glEnable(GL_POINT_SMOOTH);
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
	if (!mesh.isNull() && mesh->n_faces())
	{
		Vector3 sum(0, 0, 0);
		for (auto pi : mesh->vertices())
			sum += mesh->get_vertex_property<Vector3>("v:point")[pi];
		Vector3 c = sum / mesh->n_vertices();

		double delta = -0.5;
		c += Vector3(delta, delta, 0);

		auto points = mesh->vertex_property<Eigen::Vector3d>("v:point");
		auto colors = mesh->vertex_property<Eigen::Vector3d>("v:color");

		glBegin(GL_TRIANGLES);
		for (auto f : mesh->faces())
		{
			for (auto v : surface_mesh::Mesh::Vertex_around_face_circulator(mesh.data(), f))
			{
				auto p = points[v];
				p -= c;

				glColor3dv(colors[v].data());
				glVertex2d(p[0], p[1]);
			}
		}
		glEnd();

		// Draw boundary
		glLineWidth(2);
		glColor3d(0.5, 0.5, 0.5);
		glBegin(GL_LINE_LOOP);
		surface_mesh::Surface_mesh::Halfedge h_b, h_next;
		for (auto he : mesh->halfedges()){ if (mesh->is_boundary(he)){ h_b = he; break; } }
		h_next = mesh->next_halfedge(h_b);
		while (h_b != h_next){
			auto p = points[mesh->from_vertex(h_next)];
			p -= c;
			glVertex2d(p[0], p[1]);
			h_next = mesh->next_halfedge(h_next);
		}
		glEnd();

		glColor3d(1, 0, 0);
		glBegin(GL_POINTS);
		for (auto idx : pins)
		{
			auto p = points[surface_mesh::Surface_mesh::Vertex(idx)];
			p -= c;
			glVertex2d(p[0], p[1]);
		}
		glEnd();
	}
}

void Viewer::visualizeScalar(const string property)
{
	auto colors = mesh->vertex_property<Eigen::Vector3d>("v:color");
	auto p = mesh->vertex_property<Scalar>(property);
	for (auto v : mesh->vertices()){
		auto color = qtJetColor(p[v], 0, 1);
		colors[v] = Vector3(color.redF(), color.greenF(), color.blueF());
	}
}

void Viewer::normalizeVertexProperty(const std::string property)
{
	auto vprop = mesh->vertex_property<double>(property);
	auto V = vprop.vector();
	double minv = *std::min_element(V.begin(), V.end());
	double maxv = *std::max_element(V.begin(), V.end());
	double rangeValue = maxv - minv;
	for (auto v : mesh->vertices()) vprop[v] = (vprop[v] - minv) / rangeValue;
}

template<class Type>
Surface_mesh::Vertex_property<Type> smoothVertexProperty(surface_mesh::Surface_mesh * mesh, const std::string property, int iterations = 1, Type defaultValue = Type())
{
	Surface_mesh::Vertex_property<Type> vprop = mesh->vertex_property<Type>(property, defaultValue);

	for (int i = 0; i < iterations; i++)
	{
		std::vector<Type> newValues(mesh->n_vertices(), defaultValue);

		// average the values of neighbors
		for (auto v : mesh->vertices()){
			for (auto vj : mesh->vertices(v))
				newValues[v.idx()] += vprop[vj];
			newValues[v.idx()] /= mesh->valence(v);
		}

		// copy results back to property
		for (auto v : mesh->vertices())
			vprop[v] = newValues[v.idx()];
	}

	return vprop;
}

void Viewer::keyPressEvent(QKeyEvent *e)
{
	if (e->key() == Qt::Key_S){
		this->grabFramebuffer().save("captured.png");
	}

	if (e->key() == Qt::Key_R)
	{
		auto points = mesh->vertex_property<Eigen::Vector3d>("v:point");
		auto d = mesh->vertex_property<Scalar>("v:dist_zero");
		for (auto v : mesh->vertices())	points[v][2] = d[v];
		mesh->write("mesh.obj");

		for (auto v : mesh->vertices())	points[v][2] = 0.0;
	}

	if (e->key() == Qt::Key_D)
	{
		if (mesh.isNull()) return;

		// Get start point
		auto points = mesh->vertex_property<Eigen::Vector3d>("v:point");
		double maxDist = -DBL_MAX;
		Surface_mesh::Vertex start_v(0);
		for (auto v : mesh->vertices()){
			if (points[v].norm() > maxDist){
				maxDist = points[v].norm();
				start_v = v;
			}
		}

		QSet<Surface_mesh::Vertex> sources;
		sources << Surface_mesh::Vertex(start_v);
		auto d = GeoHeatHelper(mesh.data(), 10).getUniformDistance(sources, "v:dist_zero");
		auto bdist = mesh->vertex_property<Scalar>("v:dist_border");

		for (auto v : mesh->vertices())
		{
			double weight = d[v] + (1.0 - bdist[v]);
			d[v] = pow(weight,2);
		}

		//smoothVertexProperty<double>(mesh.data(), "v:dist_zero", 2);
		normalizeVertexProperty("v:dist_zero");
		visualizeScalar("v:dist_zero");

		update();
	}

	if (e->key() == Qt::Key_Space)
	{
		if (mesh.isNull())
		{
			trianglelib::BasicMesh<double> m;
			auto options = "pa0.0004";
			LetterGeometry letter('w', 140, 0.008);
			trianglelib::triangulatePolygon<double>(letter.points, letter.holes, m, options);

			mesh = QSharedPointer<surface_mesh::Mesh>(new surface_mesh::Mesh());

			for (auto v : m.vertices) mesh->add_vertex(Eigen::Vector3d(v[0], v[1], 0));
			for (auto t : m.triangles) mesh->add_triangle(t[0], t[1], t[2]);

			mesh->add_vertex_property<Vector3>("v:color", Vector3(1, 0, 0));

			// Pre-compute distance to border
			QSet<Surface_mesh::Vertex> sources;
			for (auto v : mesh->vertices()) if (mesh->is_boundary(v)) sources << v;
			auto dists = GeoHeatHelper(mesh.data(), 0.1).getUniformDistance(sources, "v:dist_border");
			normalizeVertexProperty("v:dist_border");
			visualizeScalar("v:dist_border");

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

				/*pins.push_back(0);
				pins.push_back(mesh->n_vertices() * 0.5);
				pins.push_back(mesh->n_vertices() * 0.25);
				pins.push_back(mesh->n_vertices() * 0.75);
				pins.push_back(mesh->n_vertices() * 0.4);
				*/

				int numPins = 5;

				GeoHeatHelper geo(mesh.data());
				auto src_points = mesh->vertex_property<bool>("v:geo_src_points", false);
				auto geo_dist = mesh->vertex_property<Scalar>("v:uniformDistance", 0);
				src_points[Vertex(0)] = true;
				for (int itr = 0; itr < numPins; itr++)
				{
					QSet<Vertex> src;
					for (auto v : mesh->vertices()) if (src_points[v]) src << v;
					geo_dist = geo.getUniformDistance(src);

					// Find furthest point, add it to source set
					double maxDist = -DBL_MAX;
					int maxIDX = -1;
					for(auto v : mesh->vertices())
					{
						if (geo_dist[v] > maxDist)
						{
							maxDist = geo_dist[v];
							maxIDX = v.idx();
						}
					}

					src_points[Vertex(maxIDX)] = true;

					pins.push_back(maxIDX);
				}

				{
					double close_weight = 1.0;

					for (auto pinned : pins)
					{
						std::vector<int> id_vector;
						id_vector.push_back(pinned);

						// Inefficient neighbors...
						QSet<int> u;
						for (auto vj : mesh->vertices(Vertex(pinned)))
						{
							u << vj.idx();
							for (auto vjj : mesh->vertices(Vertex(vj.idx())))
							{
								u << vjj.idx();
								for (auto vjjj : mesh->vertices(Vertex(vjj.idx())))
								{
									for (auto vjjjj : mesh->vertices(Vertex(vjjj.idx())))
									{
										u << vjjjj.idx();
									}
								}
							}
						}

						for (auto vj : u) id_vector.push_back(vj);

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
				for (int i = 0; i < p.cols(); i++) p.col(i) = final_points.col(i);
			}
		}

		timer = new QTimer(this);
		connect(timer, &QTimer::timeout, [&]() {
			size_t nb_points = this->mesh->n_vertices();
			auto points = mesh->vertex_property<Eigen::Vector3d>("v:point");
			Eigen::Map<ShapeOp::Matrix3X> p(points.vector().front().data(), 3, nb_points);
			int dimensions = 2;

			solver->solve(4, dimensions);

			auto final_points = solver->getPoints();
			for (int i = 0; i < p.cols(); i++) p.col(i) = final_points.col(i);

			// Smooth out spikes
			for (int i = 0; i < pins.size(); i++)
			{
				surface_mesh::Surface_mesh::Vertex vi(pins[i]);
				Vector3 pnt(0, 0, 0);
				bool isOnBoundary = mesh->is_boundary(vi);
				if (isOnBoundary)
				{
					for (auto vj : mesh->vertices(vi)) if (mesh->is_boundary(vj)) pnt += p.col(vj.idx());
					pnt /= 2;
				}
				else
				{
					for (auto vj : mesh->vertices(vi)) pnt += p.col(vj.idx());
					pnt /= mesh->valence(vi);
				}

				p.col(vi.idx()) = pnt;
			}

			update();

			for (int i = 0; i < pins.size(); i++)
			{
				if ((time.elapsed() / 1000) % 1 > 0) continue;

				if ((double(rand()) / RAND_MAX) > 0.1) continue;

				// Modify constraints
				auto c = std::dynamic_pointer_cast <ShapeOp::ClosenessConstraint>(solver->getConstraint(pins_const_ids[i]));

				Vector3 oldPos = p.col(pins[i]);
				Vector3 newPos = oldPos + (Vector3::Random() * 0.5 * (double(rand()) / RAND_MAX));
				c->setPosition(newPos);
			}


			//if (time.elapsed() > 20000) timer->deleteLater();
		});

		time.restart();
		timer->start(30);
	}
}
