#pragma once

#include <string>
#include <vector>
#include <cmath>

#include <Eigen/Geometry>

class Mesh
{
public:
    Mesh(const std::string& filename);
	Mesh(const Mesh& other);
	std::string filename;

    void loadOBJ();
	void center();

	// Vertex operations
	struct Vertex : public Eigen::Vector3d{
		Vertex(const Eigen::Vector3d & coordinate, int index = 0) : Eigen::Vector3d(coordinate), kappa(0), index(index){}
		double kappa;
		int index;

		Vertex *forward, *backward;

		double curvature(){
			double angle = phi();
			double dl = dualLength();
			return angle / dl;
		}

		double dualLength(){
			return 0.5 * ((*forward - *this).norm() + (*backward - *this).norm());
		}

		double phi(){
			auto p_next = *forward;
			auto p_prev = *backward;
			Eigen::Vector3d v = p_next - *this;
			Eigen::Vector3d u = *this - p_prev;
			return atan2(u.x() * v.y() - v.x() * u.y(), u.dot(v));
		}
	};

    std::vector<Vertex> vertices;
    int num_points;

	// Edge operations
	struct Edge{ Eigen::Vector3d tangent; };
	std::vector<Edge> edges;	
	double length(Vertex * pi, Vertex * pj){return (*pi - *pj).norm();}
	double theta(Vertex a, Vertex b){Eigen::Vector3d realtang = a - b;return atan2(realtang.y(), realtang.x());}
};
