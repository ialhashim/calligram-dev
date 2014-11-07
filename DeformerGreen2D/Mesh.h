#pragma once

#include <QString>
#include <Eigen/Geometry>

typedef double Scalar;
typedef Eigen::Vector3d Point;
typedef Eigen::Vector3d Normal;
typedef Eigen::Vector3d Vector3;
namespace{
    Eigen::Vector3d cross(Eigen::Vector3d a, Eigen::Vector3d b){ return a.cross(b); }
    double dot(Eigen::Vector3d a, Eigen::Vector3d b){ return a.dot(b); }
    double norm(Eigen::Vector3d a) { return a.norm(); }
}
#include "Surface_mesh.h"

namespace surface_mesh{
    class Mesh : public Surface_mesh{
    public:
        Mesh(QString filename, double scale = 1.0, const Point& delta = Point(0,0,0)){
            this->read_obj(filename.toStdString());

            if(scale != 1.0){
                auto & points = vertex_property<Point>("v:point");
                for(auto v : vertices()){
                    points[v] *= scale;
                    points[v] += delta;
                }
            }
        }
    };
}
