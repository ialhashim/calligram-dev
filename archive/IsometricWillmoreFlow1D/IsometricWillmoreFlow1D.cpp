#pragma warning(disable:4789)

#include <iostream>

#include "IsometricWillmoreFlow1D.h"

using namespace std;

IsometricWillmoreFlow1D::IsometricWillmoreFlow1D(Mesh * _mesh) : mesh(_mesh)
{
}

void IsometricWillmoreFlow1D::integrate(double dt)
{
    getCurvature();
    buildMassMatrix();
    buildConstraints();
    orthogonalizeConstraints();
    computeFlowDirection();
    enforceConstraints();
    integrateFlow( dt );
    recoverTangents();
    recoverPositions();
}

void IsometricWillmoreFlow1D::getCurvature(void)
{
    for (int i = 0; i < (int)mesh->vertices.size(); i++)
    {
        Mesh::Vertex & v = mesh->vertices[i];
        v.kappa = v.curvature();
    }
}

void IsometricWillmoreFlow1D::buildMassMatrix(void)
{
    int nV = (int)mesh->vertices.size();

    B = Eigen::SparseMatrix<double>(nV, nV);

    for (int i = 0; i < (int)mesh->vertices.size(); i++)
    {
        B.coeffRef(i, i) = mesh->vertices[i].dualLength();
	}
}

void IsometricWillmoreFlow1D::buildConstraints(void)
{
    int nV = (int)mesh->vertices.size();
    c.resize(3, Eigen::MatrixXd(nV, 1));

    for (int i = 0; i < (int)mesh->vertices.size(); i++)
    {
        auto & v = mesh->vertices[i];
		c[0](i, 0) = 1;
        c[1](i, 0) = v.x();
        c[2](i, 0) = v.y();
	}
}

void IsometricWillmoreFlow1D::orthogonalizeConstraints(void)
{
	double c0Norm = computeMatrixNorm(B, c[0]);
	c[0] /= c0Norm;

	// c1 = c1 - <c1, c0>c0
	double c1c0 = computeMatrixDot(B, c[1], c[0]);

	c[1] -= c1c0 * c[0];
	double c1Norm = computeMatrixNorm(B, c[1]);
	c[1] /= c1Norm;

	// c2 = c2 - <c2, c0>c0 - <c2,c1>c1
	double c2c0 = computeMatrixDot(B, c[2], c[0]);
	c[2] -= c2c0 * c[0];

	double c2c1 = computeMatrixDot(B, c[2], c[1]);
	c[2] -= c2c1 * c[1];

	double c2Norm = computeMatrixNorm(B, c[2]);
    c[2] /= c2Norm;
}

void IsometricWillmoreFlow1D::computeFlowDirection(void)
{
	size_t n = mesh->vertices.size();
    kappaDot = Eigen::MatrixXd(n, 1);

    for (int i = 0; i < (int)mesh->vertices.size(); i++)
    {
        kappaDot(i) = -2.0 * mesh->vertices[i].kappa;
    }
}

void IsometricWillmoreFlow1D::enforceConstraints(void)
{
	double coef0 = computeMatrixDot(B, kappaDot, c[0]);
	double coef1 = computeMatrixDot(B, kappaDot, c[1]);
	double coef2 = computeMatrixDot(B, kappaDot, c[2]);
	kappaDot -= coef0 * c[0];
	kappaDot -= coef1 * c[1];
	kappaDot -= coef2 * c[2];
}

void IsometricWillmoreFlow1D::integrateFlow(double dt)
{
    for (int i = 0; i < (int)mesh->vertices.size(); i++)
    {
        mesh->vertices[i].kappa += dt * kappaDot(i);
    }
}

void IsometricWillmoreFlow1D::recoverTangents(void)
{
    auto theta_i = mesh->theta(mesh->vertices.front(), mesh->vertices.back());

	for (int i = 0; i < (int)mesh->vertices.size(); i++)
	{
		auto & v = mesh->vertices[i];

        double L_i = mesh->vertices[i].dualLength();

		theta_i += L_i * v.kappa;

        mesh->edges[i].tangent = L_i * Eigen::Vector3d(cos(theta_i), sin(theta_i), 0.0);
    }
}

void IsometricWillmoreFlow1D::recoverPositions(void)
{
    Mesh::Vertex pos_i = mesh->vertices.back();

    for (int idx = 0; idx <= (int)mesh->vertices.size(); idx++)
	{
        int i = idx % mesh->vertices.size();
        int j = i - 1; if (j < 0) j += (int)mesh->edges.size();

        Eigen::Vector3d tangent = mesh->edges[j].tangent;
		pos_i += tangent;

		// Only coordinate changes
		mesh->vertices[i][0] = pos_i[0];
		mesh->vertices[i][1] = pos_i[1];
        mesh->vertices[i][2] = pos_i[2];
    }
}
