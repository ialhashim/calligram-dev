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
    void triangulatePolygon( std::vector< std::vector< Scalar > > & points, std::vector< std::vector< Scalar > > & holes, BasicMesh<Scalar>& mesh, const char * tri_switches = nullptr)
	{
        triangulateio in, out;

        // points
        in.numberofpoints = (int)points.size();
		std::vector< Scalar > localPoints( in.numberofpoints * 2, 0 );
		in.pointlist = &localPoints[0];

		for (int i = 0; i < in.numberofpoints; ++i ){
			in.pointlist[2*i + 0] = points[i][0];
			in.pointlist[2*i + 1] = points[i][1];
		}

        // segments
        in.numberofsegments = in.numberofpoints;
		std::vector< int > segmentlist( in.numberofsegments * 2, 0 );
		in.segmentlist = &segmentlist[0];

        for(int i = 0; i < in.numberofsegments; i++){
			int j = i + 1;
            if(j > in.numberofsegments - 1) j = 0;

			in.segmentlist[2*i + 0] = i;
			in.segmentlist[2*i + 1] = j;
        }

		// Holes
		in.numberofholes = 1;
		std::vector< Scalar > localHoles(in.numberofholes * 2, 0);
		in.holelist = &localHoles[0];
        if(holes.size())
        {
            int offset = (int)localPoints.size();
            localPoints.resize(offset + (2 * holes.size()), 0);

            for (int i = 0; i < holes.size(); ++i ){
				localPoints[2 * i + 0 + offset] = holes[i][0];
				localPoints[2 * i + 1 + offset] = holes[i][1];
            }

			int voffset = (int)points.size();
            segmentlist.resize(localPoints.size(), -1);
            for(int i = 0; i < holes.size(); i++){

				int idx = i + voffset;
				int j = idx + 1;

				if (i + 1 >= holes.size()) 
					j = voffset;

				segmentlist[2 * i + 0 + offset] = idx;
				segmentlist[2 * i + 1 + offset] = j;
            }

			// hole
            std::vector<double> midPoint(2,0);

            for(auto p : holes){
                midPoint[0] += p[0];
                midPoint[1] += p[1];
            }

			midPoint[0] /= holes.size();
			midPoint[1] /= holes.size();

            for (int i = 0; i < in.numberofholes; ++i ){
                in.holelist[2*i + 0] = midPoint[0];
                in.holelist[2*i + 1] = midPoint[1];
			}

			// Correct pointers
			in.numberofpoints = localPoints.size() * 0.5;
			in.pointlist = &localPoints[0];

			in.numberofsegments = segmentlist.size() * 0.5;
			in.segmentlist = &segmentlist[0];
        }

		auto defaultOpts = "pq32.5";

		triangulate(tri_switches ? tri_switches : const_cast<char*>(defaultOpts), &in, &out, NULL);

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
