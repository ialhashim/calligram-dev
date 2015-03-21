#include "Mesh.h"

#include <QTextStream>
#include <QFileInfo>
#include <QPolygonF>
#include <QPainterPath>

QPolygonF resamplePolygon(QPolygonF points, int count = 100)
{
    QPainterPath path;
    path.addPolygon(points);
    auto pathLen = path.length();
    auto stepSize = pathLen / count;

    QPolygonF newPoints;
    for(int i = 0; i < count; i++)
    {
        newPoints << path.pointAtPercent(path.percentAtLength( stepSize * i ));
    }
    return newPoints;
}

Mesh::Mesh(const std::string &filename, int resampling) : filename(filename), num_points(0)
{
    loadOBJ();

    if(resampling)
    {
        // Laplacian smoothing
        /*for(int i = 0; i < resampling; i++)
        {
            std::vector<Eigen::Vector3d> tempVertices, resultVertices(vertices.size());

            for(auto & v : vertices) tempVertices.push_back(v);

            for (int j = 0; j < (int)vertices.size(); j++)
            {
                int i = j - 1; if(i < 0) i += (int)vertices.size();
                int k = (i + 1) % vertices.size();

                Eigen::Vector3d midPoint = 0.5 * (Eigen::Vector3d(tempVertices[i]) + Eigen::Vector3d(tempVertices[k]));
                resultVertices[j] = midPoint;
            }

            for (int j = 0; j < (int)vertices.size(); j++)
            {
                vertices[j][0] = resultVertices[j][0];
                vertices[j][1] = resultVertices[j][1];
                vertices[j][2] = resultVertices[j][2];
            }
        }*/

        QPolygonF poly;
        for (int i = 0; i < (int)vertices.size(); i++)
            poly << QPointF(vertices[i][0], vertices[i][1]);
        poly << QPointF(vertices[0][0], vertices[0][1]);
        poly = resamplePolygon(poly, std::max((int)vertices.size(),resampling));

        vertices.clear();
        for (int i = 0; i < (int)poly.size(); i++)
        {
            vertices.push_back(Vertex(Eigen::Vector3d(poly[i].x(), poly[i].y(),0), i));
        }
    }

    postProcess();

    flow_vertices.push_back(vertices);
}

void Mesh::loadOBJ()
{
	auto fname = QString::fromStdString(filename);
    if(!QFileInfo(fname).exists()) return;

	QFile file(fname);
    file.open(QFile::ReadOnly | QFile::Text);
    QTextStream in(&file);
    QStringList lines = in.readAll().split('\n');
    for( auto line : lines )
    {
        auto lineTokens = line.split(" ");
        if(lineTokens.size() != 4 || lineTokens.front() != "v") continue;

		double x = lineTokens[1].toFloat(), y = lineTokens[2].toFloat(), z = lineTokens[3].toFloat();

        vertices.push_back( Vertex(Eigen::Vector3d(x,y,z), (int)vertices.size()) );
    }
}

void Mesh::postProcess()
{
    // Normalize:
    Eigen::Vector3d center(0, 0, 0);
    double rMax = 0;
    for (auto & v : vertices) center += v;
    center /= vertices.size();
    for (auto & v : vertices) v -= center;
    for (auto & v : vertices) rMax = std::max(rMax, v.norm());
    for (auto & v : vertices) v /= rMax;

    // Build connectivity:
    edges.resize(vertices.size());
    for (int i = 0; i < (int)vertices.size(); i++)
    {
        int j = (i + 1) % (vertices.size());
        vertices[i].forward = &vertices[j];

        int k = i - 1; if (k < 0) k += (int)vertices.size();
        vertices[i].backward = &vertices[k];
    }
}

void Mesh::center()
{
	Eigen::Vector3d c(0, 0, 0);
	double L = 0;

	for (int i = 0; i < (int)vertices.size(); i++)
	{
		Eigen::Vector3d a = vertices[i];
		Eigen::Vector3d b = vertices[(i+1)%vertices.size()];

		Eigen::Vector3d midpoint = (a + b) / 2.;
		double length = (a - b).norm();

		c += length * midpoint;
		L += length;
	}

	c /= L;

	for (int i = 0; i < (int)vertices.size(); i++)
	{
		auto & v = vertices[i];
		v -= c;
	}
}
