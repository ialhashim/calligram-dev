#include "Mesh.h"

#include <QTextStream>
#include <QFileInfo>

Mesh::Mesh(const std::string &filename) : filename(filename), num_points(0)
{
    loadOBJ();
}

Mesh::Mesh(const Mesh &other) : filename(other.filename), num_points(0)
{
    loadOBJ();
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
