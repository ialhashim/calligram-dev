#pragma once

#include "Mesh.h"
#include <Eigen/Sparse>

class IsometricWillmoreFlow1D
{
public:
	IsometricWillmoreFlow1D(Mesh * mesh);
	Mesh * mesh;

    void integrate(double dt = 1e-2);

protected:
    void getCurvature(void);
    void buildMassMatrix(void);
    void buildConstraints(void);
    void orthogonalizeConstraints(void);
    void computeFlowDirection(void);
    void enforceConstraints(void);
    void integrateFlow(double dt);
    void recoverTangents(void);
    void recoverPositions(void);

    Eigen::SparseMatrix<double> B;      // diagonal mass matrix on primal 0-forms
    std::vector< Eigen::MatrixXd > c;   // constraint vectors
    Eigen::MatrixXd kappaDot;           // flow direction

    double computeMatrixDot(const Eigen::SparseMatrix<double>& M, const Eigen::MatrixXd& v1, const Eigen::MatrixXd& v2){
        auto Mv = M * v2;
        auto res = v1.transpose() * Mv;
        return res(0, 0);
    }

    double computeMatrixNorm(const Eigen::SparseMatrix<double>& M, const Eigen::MatrixXd& v){
        return sqrt(computeMatrixDot(M, v, v));
    }
};

