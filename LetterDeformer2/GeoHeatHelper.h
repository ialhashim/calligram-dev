#pragma once
#include <float.h>

// Eigen matrix library
#include <Eigen/Core>
#include <Eigen/Sparse>

using namespace std;
using namespace Eigen;
using namespace surface_mesh; 

#define qRanged(min, v, max) ( qMax(min, qMin(v, max)) )
#define EPSILON 1e-12
typedef Eigen::SimplicialLLT< SparseMatrix<double> > GeoHeatSolver;

typedef Surface_mesh::Vertex Vertex;

static uint qHash(const Vertex &key){ return qHash(key.idx()); }

class GeoHeatHelper{

public:
	typedef Surface_mesh::Vertex_property<Vector3>  Vector3VertexProperty;  ///< An Vector3 associated to a vertex.
	typedef Surface_mesh::Vertex_property<Scalar>   ScalarVertexProperty;   ///< A scalar associated to a vertex.
	typedef Surface_mesh::Edge_property<Scalar>     ScalarEdgeProperty;     ///< A scalar associated to an edge.
	typedef Surface_mesh::Face_property<Vector3>    Vector3FaceProperty;    ///< A Vector3 associated to a face.
	typedef Surface_mesh::Face_property<Scalar>     ScalarFaceProperty;     ///< A scalar associated to a face.
	
	Vector3VertexProperty points;
	ScalarVertexProperty  varea;
	Vector3VertexProperty vnormal;
	Vector3FaceProperty   fnormal; 
	ScalarFaceProperty    farea; 
	ScalarEdgeProperty    elenght;  

    ScalarVertexProperty vfunction;
    ScalarVertexProperty vcot;
    ScalarVertexProperty u;
    ScalarVertexProperty vdiv;
    ScalarEdgeProperty ecot;
    Vector3FaceProperty fgradient;

    SparseMatrix<Scalar>    Lc_;
	GeoHeatSolver heat_flow, poisson_solver;

	Scalar t_factor;

	Surface_mesh* mesh;

private:
	ScalarEdgeProperty computeEdgeLengths(std::string property = "e:length"){
		elenght = mesh->edge_property<Scalar>(property, 0.0f);
		for(auto eit : mesh->edges())
			elenght[eit] = mesh->edge_length(eit);
		return elenght;
	}
	ScalarVertexProperty computeVertexVoronoiArea(const std::string property = "v:varea"){
		varea = mesh->vertex_property<Scalar>(property);
		Scalar a;
		Vertex v0, v1, v2;
		for(auto v : mesh->vertices())
			varea[v] = 0.0;
		for(auto f : mesh->faces()){
			Surface_mesh::Vertex_around_face_circulator vfit = mesh->vertices(f);
			v0 = *vfit;
			v1 = *++vfit;
			v2 = *++vfit;

			// compute area
			a = 0.5 * (points[v1] - points[v0]).cross(points[v2] - points[v0]).norm();

			// distribute area to vertices
			varea[v0] += a / 3.0;
			varea[v1] += a / 3.0;
			varea[v2] += a / 3.0;
		}
		return varea;
	}

	Surface_mesh::Face_property<Scalar> computeFaceAreas(std::string property = "f:area"){
		farea = mesh->face_property<Scalar>(property);

		Surface_mesh::Face_iterator fit, fend = mesh->faces_end();
		Surface_mesh::Vertex_around_face_circulator vit, vend;
		QVector<Vector3> pnts;

		for (fit = mesh->faces_begin(); fit != fend; ++fit){
			// Collect points of face
			pnts.clear(); vit = vend = mesh->vertices(*fit);
			do{ pnts.push_back(points[*vit]); } while (++vit != vend);

			farea[*fit] = 0.5 * (pnts[1] - pnts[0]).cross(pnts[2] - pnts[0]).norm();
		}

		return farea;
	}


public:
	GeoHeatHelper(Surface_mesh* mesh, double t_factor = 1.0) : mesh(mesh), t_factor(t_factor)
	{
		// Prepare mesh properties
		points  = mesh->vertex_property<Point>("v:point");
        varea   = this->computeVertexVoronoiArea();
        farea   = this->computeFaceAreas();
		elenght = this->computeEdgeLengths();
		mesh->update_face_normals();
		fnormal = mesh->face_property<Vector3>("f:normal");

        // Heat method properties
        u           = mesh->vertex_property<Scalar> ("v:heat", 0);
        vdiv        = mesh->vertex_property<Scalar> ("v:divergence", 0);
        vcot        = mesh->vertex_property<Scalar> ("v:cotan", 0);
        vfunction   = mesh->vertex_property<Scalar> ("v:function", 0);
        ecot        = mesh->edge_property<Scalar>   ("e:cotan", 0);
        fgradient   = mesh->face_property<Vector3>  ("f:gradient", Vector3(0,0,0));
    }

    ~GeoHeatHelper()
    {
        cleanUp( true );
    }

    void precompute()
    {
        Lc_ = Lc();

		heat_flow.compute(A() + (t() * Lc_));
        poisson_solver.compute( Lc_ );
    }

	ScalarVertexProperty getUniformDistance(const QSet<Vertex> &source, const string property = "v:uniformDistance")
    {
		if(source.empty()) return ScalarVertexProperty();

        if(Lc_.size() == 0) precompute();

        // 0) Set source vertices
        VectorXd u0_ = u0( source );

        // 1) Compute heat flow for time 't'
        set( heat_flow.solve( u0_ ), "v:heat" );

        // 2) Evaluate vector field X
        gradientFaces();

        // 3) Compute distance function (solve Poisson equation)
        VectorXd d = divergenceVertices();
        set( poisson_solver.solve( d ), "v:heat_distance");

		return unifromDistance(property);
    }

    Scalar t()
    {
        // t = (avg edge length) ^ 2
        double avg = 0.0;
        int count = 0;

        for(auto e: mesh->edges()){
            if(elenght[e] > 0.0){
                avg += elenght[e];
                count++;
            }
        }

        avg /= count;
        return (avg * avg) * t_factor;
    }

    VectorXd u0(const QSet<Vertex> & vidx)
    {
        mesh->remove_vertex_property(vfunction);
        vfunction = mesh->vertex_property<Scalar> ("v:function", 0);

        // u0 on source vertices = 1.0 and zero elsewhere
        VectorXd u_0 = VectorXd::Zero(mesh->n_vertices());
        for(auto v : vidx)
            u_0[v.idx()] = vfunction[v] = 1.0;

        return u_0;
    }

    SparseMatrix<Scalar> A()
    {
        // A is diagonal matrix of vertex areas
        typedef Eigen::Triplet<double> T;
        std::vector< T > A_elements;

        for(auto v : mesh->vertices())
            A_elements.push_back( T(v.idx(), v.idx(), qMax(varea[v], EPSILON) ) );

        SparseMatrix<Scalar> matA(mesh->n_vertices(), mesh->n_vertices());
        matA.setFromTriplets(A_elements.begin(), A_elements.end());
        return matA;
    }

    SparseMatrix<Scalar> Lc()
    {
        // Efficient sparse matrix construction
        typedef Eigen::Triplet<double> T;
        std::vector< T > L_c;

        // Fill as cotan operator
        for(auto e : mesh->edges())
		{
            Scalar cot_alpha = 0, cot_beta = 0;

            Vertex vi = mesh->vertex(e, 0);
            Vertex vj = mesh->vertex(e, 1);

            Vertex v_a = mesh->to_vertex(mesh->next_halfedge(mesh->halfedge(e, 0)));
            if(has_halfedge(v_a, vj)) cot_alpha = (points[vi]-points[v_a]).dot(points[vj]-points[v_a]) / (points[vi]-points[v_a]).cross(points[vj]-points[v_a]).norm();

            Vertex v_b = mesh->to_vertex(mesh->next_halfedge(mesh->halfedge(e, 1)));
            if(has_halfedge(v_b, vi)) cot_beta  = (points[vi]-points[v_b]).dot(points[vj]-points[v_b]) / (points[vi]-points[v_b]).cross(points[vj]-points[v_b]).norm();

            Scalar cots = (0.5 * (cot_alpha + cot_beta));

            if(abs(cots) == 0) continue;

            L_c.push_back(T(vi.idx(), vj.idx(), -cots));
            L_c.push_back(T(vj.idx(), vi.idx(), -cots));
            L_c.push_back(T(vi.idx(), vi.idx(), cots + EPSILON));
            L_c.push_back(T(vj.idx(), vj.idx(), cots + EPSILON));

            // Just for record
            vcot[vi] += cots;
            vcot[vj] += cots;
            ecot[e] = cots;
        }

        // Initialize a sparse matrix
        SparseMatrix<Scalar> Lc_mat(mesh->n_vertices(), mesh->n_vertices());
        Lc_mat.setFromTriplets(L_c.begin(), L_c.end());
        return Lc_mat;
    }

    void gradientFaces(bool isNormalizeNegateGradient = true)
    {
        // Compute gradient on faces
        for(auto f : mesh->faces())
		{
            Vector3 vsum(0,0,0);

            Surface_mesh::Halfedge_around_face_circulator h(mesh, f), hend = h;
            do{
                Vector3 ei = points[mesh->from_vertex(*h)] - points[mesh->from_vertex(mesh->prev_halfedge(*h))];
                Vertex i = mesh->to_vertex(*h);
                vsum += u[i] * fnormal[f].cross(ei);
            }while (++h != hend);

            fgradient[f] = ( 1.0 / (2.0 * farea[f]) ) * vsum;

            if(isNormalizeNegateGradient)
                fgradient[f] = (-fgradient[f]).normalized();
        }
    }

    VectorXd divergenceVertices()
    {
        VectorXd d(mesh->n_vertices());

        for(auto i : mesh->vertices())
		{
            double sum_j = 0.0;

            Surface_mesh::Halfedge_around_vertex_circulator j(mesh, i), hend = j;
            do {
                if(!mesh->is_valid(mesh->face(*j))) continue; // todo: handle this?

                // Face gradient
                Vector3 Xj = fgradient[ mesh->face(*j) ];

                // Vertex position
                Vector3 pi = points[mesh->from_vertex(*j)];
                Vector3 p1 = points[mesh->to_vertex(*j)];
                Vector3 p2 = points[mesh->to_vertex(mesh->next_halfedge(*j))];

                // Incident edges
                Vector3 e1 = p1 - pi;
                Vector3 e2 = p2 - pi;

                // Angles
                double theta1 = acos( ((p1-p2).normalized()).dot((pi-p2).normalized()) );
                double theta2 = acos( ((p2-p1).normalized()).dot((pi-p1).normalized()) );
                double cot1 = 1.0 / tan(theta1);
                double cot2 = 1.0 / tan(theta2);

                sum_j += (cot1 * (e1).dot(Xj)) + (cot2 * (e2).dot(Xj));

            } while(++j != hend);

            vdiv[i] = 0.5 * sum_j;

            d(i.idx()) = vdiv[i];
        }

        return d;
    }

    VectorXd toEigenVector(const ScalarVertexProperty & vproperty){
        VectorXd V( mesh->n_vertices() );
        for(auto i : mesh->vertices()) V(i.idx()) = vproperty[i];
        return V;
    }

    void set(const VectorXd & V, const string property){
        ScalarVertexProperty vprop = mesh->vertex_property<Scalar>(property);
		for (auto v : mesh->vertices())
            vprop[v] = V[v.idx()];
    }

    bool has_halfedge(Vertex start, Vertex end){
        auto h = mesh->halfedge(start);
        auto hh = h;
        if (h.is_valid()){
            do{
                if (mesh->to_vertex(h) == end) return true;
                h = mesh->cw_rotated_halfedge(h);
            } while (h != hh);
        }
        return false;
    }

    ScalarVertexProperty unifromDistance( const string property = "v:uniformDistance" ){
        ScalarVertexProperty vprop = mesh->vertex_property<Scalar>(property);
        ScalarVertexProperty dist = mesh->vertex_property<Scalar>("v:heat_distance");

        double minDist = DBL_MAX, maxDist = -DBL_MAX;
		for (auto v : mesh->vertices()){
            minDist = qMin(dist[v], minDist);
            maxDist = qMax(dist[v], maxDist);
        }
        double range = maxDist - minDist;

		for (auto v : mesh->vertices())
            vprop[v] = 1 - ((dist[v] - minDist) / range);

        return vprop;
    }

	// Assuming "unifromDistance(pname)" has been called with default name
	std::vector<Vertex> shortestVertexPath(const Vertex & toVertex)
	{
		ScalarVertexProperty dists = mesh->vertex_property<Scalar>("v:uniformDistance");

		std::vector<Vertex> path;
		path.push_back(toVertex);

		// back track
		auto h = mesh->halfedge( toVertex );
		Vertex curV = toVertex;

		while( dists[mesh->to_vertex(h)] != 0.0 ) // safe?
		{
			double minDist = DBL_MAX;
			for (auto vj : mesh->vertices(curV))
			{
				double dist = dists[ vj ];
				if(dist < minDist){
					minDist = dist;
					curV = vj;
				}
			}

			if(curV != path.back())
				path.push_back(curV);
			else
				break;
		}

		std::reverse(path.begin(), path.end());
		return path;
	}

    void cleanUp(bool isAll = false)
    {
		auto heat = mesh->vertex_property<Scalar>("v:heat_distance");
		auto udist = mesh->vertex_property<Scalar>("v:uniformDistance");
                
        ecot        = mesh->edge_property<Scalar>   ("e:cotan", 0);
        fgradient   = mesh->face_property<Vector3>  ("f:gradient", Vector3(0,0,0));

        mesh->remove_vertex_property(u);
        mesh->remove_vertex_property(vdiv);
        mesh->remove_vertex_property(vcot);
        mesh->remove_vertex_property(vfunction);
        mesh->remove_edge_property(ecot);
        mesh->remove_face_property(fgradient);
        mesh->remove_vertex_property(heat);
        
        if(isAll) mesh->remove_vertex_property(udist);
    }
};
