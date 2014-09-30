#pragma once

#include <vector>

#define REAL double
#include "triangle.h"

namespace trianglelib{
	template<typename Scalar>
	struct BasicMesh{
		std::vector< std::vector<Scalar> > vertices;
		std::vector< std::vector<int> > triangles;
	};

	template<typename Scalar>
	void triangulatePolygon( std::vector< std::vector< Scalar > > & points, BasicMesh<Scalar>& mesh, char tri_switches[] = "pq" )
	{
		triangulateio in, out, vorout;

		in.numberofpoints = (int)points.size();
		std::vector< Scalar > localPoints( in.numberofpoints * 2, 0 );
		in.pointlist = &localPoints[0];

		for (int i = 0; i < in.numberofpoints; ++i ){
			in.pointlist[2*i + 0] = points[i][0];
			in.pointlist[2*i + 1] = points[i][1];
		}

		in.numberofsegments = in.numberofpoints;
		std::vector< int > segmentlist( in.numberofsegments * 2, 0 );
		in.segmentlist = &segmentlist[0];

		for(int i = 0; i < in.numberofsegments; i++){
			int j = i + 1;
			if(j > in.numberofsegments - 1) j = 0;

			in.segmentlist[2*i + 0] = i;
			in.segmentlist[2*i + 1] = j;
		}

		triangulate(tri_switches, &in, &out, NULL);

		// Add mesh vertices
		for(int i = 0; i < out.numberofpoints; i++)
		{
			std::vector<Scalar> vertex(2,0);
			vertex[0] = out.pointlist[i * 2 + 0];
			vertex[1] = out.pointlist[i * 2 + 1];
			mesh.vertices.push_back( vertex );
		}

		// Add mesh triangles
		for(int i = 0; i < out.numberoftriangles; i++)
		{
			std::vector<int> triangle(3,0);
			triangle[0] = out.trianglelist[i * 3 + 0];
			triangle[1] = out.trianglelist[i * 3 + 1];
			triangle[2] = out.trianglelist[i * 3 + 2];

			mesh.triangles.push_back( triangle );
		}
	}
}
