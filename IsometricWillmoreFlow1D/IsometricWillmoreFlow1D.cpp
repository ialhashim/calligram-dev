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
	integrateFlow(dt);
	recoverTangents();
	recoverPositions();
}

void IsometricWillmoreFlow1D::getCurvature(void)
{
	int n = mesh->vertices.size();

	double sum_kappa = 0;

	for (int i = 0; i < (int)n; i++)
	{
		auto & v = mesh->vertices[i];
		auto kappa = v.curvature();
        v.kappa = kappa;

		sum_kappa += kappa;
    }

	//printf("sum_kappa = %f\n", sum_kappa);
}

void IsometricWillmoreFlow1D::buildMassMatrix(void)
{
    int nV = (int)mesh->vertices.size();

    B = Eigen::SparseMatrix<double>(nV, nV);

    for (int i = 0; i < (int)mesh->vertices.size(); i++)
    {
        auto & v = mesh->vertices[i];
        B.coeffRef(i, i) = v.dualLength();
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

    /*
	//check orthogonality
	DenseMatrix<double> check1 = c[1].transpose() * (B * c[0]);
	double r = check1(0, 0);
    DenseMatrix<double> check2 = c[2].transpose() * (B * c[0]);
	double r1 = check2(0, 0);
	DenseMatrix<double> check3 = c[2].transpose() * (B * c[1]);
	double r2 = check3(0, 0);
	double eps = 10e-8;
	assert(fabs(r) < eps && fabs(r1) < eps && fabs(r2) < eps);
    */
}

void IsometricWillmoreFlow1D::computeFlowDirection(void)
{
	size_t n = mesh->vertices.size();
    kappaDot = Eigen::MatrixXd(n, 1);

	double sum_kappaDot = 0;

    for (int i = 0; i < (int)mesh->vertices.size(); i++)
    {
        auto & v = mesh->vertices[i];
        kappaDot(i) = -2.0 * v.kappa;

		sum_kappaDot += kappaDot(i);
	}
	//printf("computeFlowDirection = %f\n", sum_kappaDot);
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
	double sum_kappa = 0;

    for (int i = 0; i < (int)mesh->vertices.size(); i++)
    {
        auto & v = mesh->vertices[i];
		v.kappa += dt * kappaDot(i);
		sum_kappa += v.kappa;
	}

	//printf("\nintegrateFlow = %f\n", sum_kappa);
}

void IsometricWillmoreFlow1D::recoverTangents(void)
{
	auto theta_i = mesh->theta(mesh->vertices.front(), mesh->vertices.back());

	for (int i = 0; i < (int)mesh->vertices.size(); i++)
	{
		auto & v = mesh->vertices[i];

		double L_i = v.dualLength();
		theta_i += L_i * v.kappa;

		double x = cos(theta_i);
		double y = sin(theta_i);

		mesh->edges[i].tangent = L_i * Eigen::Vector3d(x, y, 0.0);
	}

	//printf("\nrecoverTangents = %f\n", theta_i);
}

void IsometricWillmoreFlow1D::recoverPositions(void)
{
	int idx_start = mesh->vertices.size() - 1;

	Mesh::Vertex pos_start = mesh->vertices[idx_start];
	Mesh::Vertex pos_i = mesh->vertices[idx_start];

	//cout << "start idx = " << idx_start << " pos = " << pos_i.x() << "," << pos_i.y() << "\n";

	for (int i = 0; i < (int)mesh->vertices.size(); i++)
	{
		int j = i - 1;
		if (j < 0) j += mesh->edges.size();
		Eigen::Vector3d tangent = mesh->edges[j].tangent;
		pos_i += tangent;

		// Only coordinate changes
		mesh->vertices[i][0] = pos_i[0];
		mesh->vertices[i][1] = pos_i[1];
		mesh->vertices[i][2] = pos_i[2];
		//cout << "idx = " << i << " Pos = " << pos_i.x() << "," << pos_i.y() << " Tangent idx " << j << " = " << mesh->edges[i].tangent.x() << "\n";
	}

	auto pp = pos_i;
}
