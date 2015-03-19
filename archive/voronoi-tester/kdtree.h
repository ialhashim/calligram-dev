#pragma once

#include "nanoflann.hpp"

template <typename Vector3, typename T>
struct PointCloud{
    std::vector<Vector3>  pts;
    inline size_t kdtree_get_point_count() const { return pts.size(); }
    inline T kdtree_distance(const T *p1, const size_t idx_p2, size_t size) const{
        const T d0=p1[0]-pts[idx_p2].x();
        const T d1=p1[1]-pts[idx_p2].y();
        //const T d2=p1[2]-pts[idx_p2].z();
        size = size;
        return (d0*d0 + d1*d1 /*+ d2*d2*/);
    }
    inline T kdtree_get_pt(const size_t idx, int dim) const{
        if (dim==0) return pts[idx].x();
        else return pts[idx].y();
    }
    template <class BBOX>
    bool kdtree_get_bbox(BBOX &) const { return false; }
};

template<typename PointType>
struct KDTree2D{
    typedef nanoflann::KDTreeSingleIndexAdaptor< nanoflann::L2_Simple_Adaptor<double, PointCloud<PointType, double> >, PointCloud<PointType, double>, 2 /* dim */ > my_kd_tree;
    my_kd_tree * tree;
    PointCloud<PointType, double> cloud;

    KDTree2D(std::vector<PointType> points){
        for(auto p : points) cloud.pts.push_back(p);
        tree = new my_kd_tree(2, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
        tree->buildIndex();
    }

    struct KDResult{ int idx; double dist; KDResult(int idx = -1, double dist = DBL_MAX):idx(idx),dist(dist){ } };
    std::vector<KDResult> search_k(PointType p, int k = 1)
    {
        std::vector<KDResult> res;

        std::vector<size_t> ret_index(k);
        std::vector<double> out_dist(k);

        double point[] = { p.x(), p.y() };
        tree->knnSearch(point, k, &ret_index[0], &out_dist[0]);

        for(size_t i = 0; i < k; i++)
            res.push_back( KDResult(ret_index[i], out_dist[i]) );

        return res;
    }
};
