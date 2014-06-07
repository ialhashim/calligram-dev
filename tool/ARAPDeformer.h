#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SVD>
#include <Eigen/Geometry>
using namespace Eigen;


#include "trianglelib.h"

typedef double Scalar;
typedef Eigen::Vector3d Point;
typedef Eigen::Vector3d Normal;
namespace{
	Eigen::Vector3d cross(Eigen::Vector3d a, Eigen::Vector3d b){ return a.cross(b); }
	double dot(Eigen::Vector3d a, Eigen::Vector3d b){ return a.dot(b); }
	double norm(Eigen::Vector3d a) { return a.norm(); }
}
#include "Surface_mesh.h"
using namespace surface_mesh;

struct ARAPDeformer{

public:

	std::vector<Matrix3d> R;
	std::vector<Vector3d> OrigMesh;
	std::vector<VectorXd> xyz;
	std::vector<VectorXd> b;

	// Frequently used
	int nVerts;
	Surface_mesh::Vertex_property<Vector3d> points;
	Surface_mesh::Vertex_property<Vector3d> normals;
	Surface_mesh::Vertex_iterator vit, vend;
	Surface_mesh::Vertex_property< std::map<Surface_mesh::Vertex, double> > wij_weight;
	Surface_mesh::Vertex_property< bool > isAnchorPoint, isControlPoint;
	Surface_mesh::Vertex_around_vertex_circulator vvit, vvend;

	SparseMatrix<double> At;
	Eigen::SimplicialLLT< SparseMatrix<double> > solver;
	bool isSolverReady;

	Surface_mesh * mesh;

	// Control points
	void SetAnchor( const Surface_mesh::Vertex & v ){
		isAnchorPoint[v] = true;
	}

	void UpdateControl( const Surface_mesh::Vertex & v, const Eigen::Vector3d & newPos ){
		isControlPoint[v] = true;
		points[v] = newPos;
		isSolverReady = false;
	}

	void SetControl( const Surface_mesh::Vertex & v ){
		isControlPoint[v] = true;
		isSolverReady = false;
	}

	void PushControl( const Surface_mesh::Vertex & v, const Eigen::Vector3d & newPos )
	{
		points[v] += newPos;
	}

	void SetControlPosition( const Surface_mesh::Vertex & v, const Eigen::Vector3d & newPos )
	{
		points[v] = newPos;
	}

	void ClearAnchors(){
		for (vit = mesh->vertices_begin(); vit != vend; ++vit)
			isAnchorPoint[vit] = false;
		isSolverReady = false;
	}

	void ClearControl(){
		for (vit = mesh->vertices_begin(); vit != vend; ++vit)
			isControlPoint[vit] = false;
		isSolverReady = false;
	}

	void ClearAll(){
		ClearAnchors();
		ClearControl();
	}

	Surface_mesh::Halfedge find_halfedge(Surface_mesh::Vertex start, Surface_mesh::Vertex end){
		Surface_mesh::Halfedge h  = mesh->halfedge(start);
		const Surface_mesh::Halfedge hh = h;
		if (h.is_valid()){
			do{
				if (mesh->to_vertex(h) == end)return h;
				h = mesh->cw_rotated_halfedge(h);
			} while (h != hh);
		}
		return Surface_mesh::Halfedge();
	}

	ARAPDeformer(Surface_mesh * usingMesh) : mesh(usingMesh), isSolverReady(false)
	{
		// Frequently used
		nVerts = mesh->n_vertices();
		vend = mesh->vertices_end();
		points = mesh->vertex_property<Vector3d>("v:point");
		normals = mesh->vertex_property<Vector3d>("v:normal");
		isAnchorPoint = mesh->vertex_property< bool >("v:anchorPoint", false);
		isControlPoint = mesh->vertex_property< bool >("v:controlPoint", false);
		wij_weight = mesh->vertex_property< std::map<Surface_mesh::Vertex, double> >("v:wij_weight");
	}

	~ARAPDeformer()
	{
		mesh->remove_vertex_property(isAnchorPoint);
		mesh->remove_vertex_property(isControlPoint);
		mesh->remove_vertex_property(wij_weight);
	}

	Vector3d cross(const Vector3d& a, const Vector3d&b){ return a.cross(b); }

	void ComputeCotWeights()
	{
		Vector3d p, p1, p2, p3;
		Surface_mesh::Vertex pid, p1id, p2id, p3id;

		double wij = 0;

		for (vit = mesh->vertices_begin(); vit != vend; ++vit)
		{
			// find alpha and beta (LSO paper), then calculate cotangent weights
			// alpha is P_P1_P2, beta is P_P3_P2
			vvit = vvend = mesh->vertices(vit); 
			do{ 
				p  = points[vit];	pid  = vit;
				p2 = points[vvit];	p2id = vvit; --vvit;
				p3 = points[vvit];	p3id = vvit; ++vvit;  ++vvit;
				p1 = points[vvit];	p1id = vvit; --vvit;

				// wij = 1/2 * (cot(alpha) + cot(beta)), for boundary edge, there is only one such edge
				wij = 0;

				if(!mesh->is_boundary(find_halfedge(pid, p2id))
					&&!mesh->is_boundary(find_halfedge(p2id, pid))) // not a boundary edge
					wij = dot(Vector3d(p-p1),Vector3d(p2-p1)) / cross(Vector3d(p-p1),Vector3d(p2-p1)).norm() + dot(Vector3d(p-p3),Vector3d(p2-p3)) / cross(Vector3d(p-p3),Vector3d(p2-p3)).norm();
				else // boundary edge, only have one such angle
				{
					if(p1id == p3id) // two angles are the same, e.g. corner of a square
					{
						wij = dot(Vector3d(p-p1),Vector3d(p2-p1)) / cross(Vector3d(p-p1),Vector3d(p2-p1)).norm();
					}
					else // find the angle not on the boundary
					{
						if(!mesh->is_boundary(find_halfedge(pid, p1id))
							&&!mesh->is_boundary(find_halfedge(p1id, pid)))
							wij = dot(Vector3d(p-p1),Vector3d(p2-p1)) / cross(Vector3d(p-p1),Vector3d(p2-p1)).norm();
						else
							wij = dot(Vector3d(p-p3),Vector3d(p2-p3)) / cross(Vector3d(p-p3),Vector3d(p2-p3)).norm();
					} 
				}

				wij_weight[vit][vvit] = wij / 2.0;

			} while(++vvit != vvend);
		}
	}

	void BuildAndFactor(){
		ComputeCotWeights();

		// Initialize
		R.clear();			R.resize(nVerts, Matrix3d::Identity());
		xyz.clear();		xyz.resize(3, VectorXd::Zero(nVerts));
		b.clear();			b.resize(3, VectorXd::Zero(nVerts));
		OrigMesh.clear();	OrigMesh.resize(nVerts, Vector3d::Zero());

		// L matrix, n by n, cotangent weights
		typedef Eigen::Triplet<Scalar> T;
		std::vector< T > L;
		for (vit = mesh->vertices_begin(); vit != vend; ++vit)
		{
			int i = Surface_mesh::Vertex(vit).idx();

			OrigMesh[i] = Vector3d(points[vit][0],points[vit][1],points[vit][2]);

			double weight = 0.0;

			if(!(isAnchorPoint[vit] || isControlPoint[vit]))
			{
				vvit = vvend = mesh->vertices(vit);	
				do{
					int j = (*vvit).idx();
					weight += wij_weight[vit][vvit];
					L.push_back(T(i, j, -wij_weight[vit][vvit]));
				} while(++vvit != vvend);
			}
			else
				weight = 1.0;

			L.push_back(T(i, i, weight));
		}

		SparseMatrix<double> A(nVerts, nVerts);
		A.setFromTriplets(L.begin(), L.end());

		At = A.transpose();

		// FACTOR:
		solver.compute(At * A);

		isSolverReady = true;
	}

	void SVDRotation()
	{
		Matrix3d eye = Matrix3d::Identity();

		for (vit = mesh->vertices_begin(); vit != vend; ++vit)
		{
			int i = Surface_mesh::Vertex(vit).idx();
			int valence = mesh->valence(vit);
			int degree = 0;

			MatrixXd P(3, valence), Q(3, valence);

			vvit = vvend = mesh->vertices(vit);
			do{
				// eij = pi - pj, pi is v_it, pj is vv_it, including weights wij
				int j = (*vvit).idx();

				P.col(degree) = (OrigMesh[i] - OrigMesh[j]) * wij_weight[vit][vvit];
				Q.col(degree++) = (	Vector3d(xyz[0][i], xyz[1][i], xyz[2][i]) -
					Vector3d(xyz[0][j], xyz[1][j], xyz[2][j]));
			} while(++vvit != vvend);

			// Compute the 3 by 3 covariance matrix:
			// actually S = (P * W * Q.t()); W is already considerred in the previous step (P=P*W)
			MatrixXd S = (P * Q.transpose());

			// Compute the singular value decomposition S = UDV.t
			JacobiSVD<MatrixXd> svd(S, ComputeThinU | ComputeThinV); // X = U * D * V.t()

			MatrixXd V = svd.matrixV();
			MatrixXd Ut = svd.matrixU().transpose();

			eye(2,2) = (V * Ut).determinant();	// remember: Eigen starts from zero index

			// V*U.t may be reflection (determinant = -1). in this case, we need to change the sign of
			// column of U corresponding to the smallest singular value (3rd column)
			R[i] = (V * eye * Ut); //Ri = (V * eye * U.t());
		}
	}

	void Deform( int ARAPIteration /*= 1*/ )
	{
		if(!isSolverReady)
			BuildAndFactor();

		// ARAP iteration
		for(int iter = 0; iter <= ARAPIteration; iter++)
		{	
			// update vector b3 = wij/2 * (Ri+Rj) * (pi - pj), where pi and pj are coordinates of the original mesh
			for (vit = mesh->vertices_begin(); vit != vend; ++vit)
			{
				int i = Surface_mesh::Vertex(vit).idx();

				Vector3d p (points[vit][0], points[vit][1], points[vit][2]);

				if(!(isAnchorPoint[vit] || isControlPoint[vit]))
				{
					p = Vector3d::Zero(); // Set to zero

					// Collect neighbors
					vvit = vvend = mesh->vertices(vit);	
					do{ 
						int j = (*vvit).idx();

						Vector3d pij = OrigMesh[i] - OrigMesh[j];
						Vector3d RijPijMat = ((R[i] + R[j]) * pij);

						p += RijPijMat * (wij_weight[vit][vvit] / 2.0);

					} while(++vvit != vvend);
				}

				// Set RHS
				for(int k = 0; k < 3; k++)
					b[k][i] = p[k];
			}

			// SOLVE for x, y, and z
			for(int k = 0; k < 3; k++)
				xyz[k] = solver.solve(At * b[k]);

			// if iter = 0, just means naive Laplacian Surface Editing (Ri is identity matrix)
			if(iter > 0) SVDRotation();
		}

		// update vertex coordinates 
		for (vit = mesh->vertices_begin(); vit != vend; ++vit)
		{
			int i = Surface_mesh::Vertex(vit).idx();
			points[vit] = Eigen::Vector3d (xyz[0][i], xyz[1][i], xyz[2][i]);
		}
	}

};
