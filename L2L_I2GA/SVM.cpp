#include "SVM.h"
#include "MatrixElement.h"
#include <iostream>
#include <vector>
#include <cmath>
#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <time.h>
#include <random>
using namespace Eigen;
GeneralizedSelfAdjointEigenSolver<MatrixXd> ges;
VectorXd D;
double E;

std::random_device Rd;
std::default_random_engine eng2(Rd());
std::uniform_real_distribution<double> distr2(0.0, 1.0);


//=============================================================================
SVM::SVM() : me()
{
	bmax = 10;
	int ndb = 100;
	Hmatrix = MatrixXd::Zero(ndb, ndb);
	Nmatrix = MatrixXd::Zero(ndb, ndb);
}
//=============================================================================
int SVM::CheckOverlap(vector<Vector3d>& Basis)
{
	int itr = Basis.size() - 1;
	float vnorm, vdotv;
	MatrixXd Nmatrix = NormMatrix(Basis);
	vnorm = Nmatrix(itr, itr);
	if (vnorm < 1e-9) return 0;
	if (itr > 0) {
		for (int i = 0; i < itr; i++) {
			vdotv = Nmatrix(itr, i) / sqrt(Nmatrix(i, i) * Nmatrix(itr, itr));
			if (vdotv > 0.9997) {
				return 0;
			}
		}
	}
	return 1;
}
//=============================================================================
double EigenValuesEquation(int itr, VectorXd D, VectorXd q, double aa, double xx)
{
	double vv = 1;
	double ww = 1;
	double yy = 0;
	double zz = 0;

	for (int n1 = 0; n1 < itr; n1++)
	{
		vv = vv * (D(n1) - xx);
		ww = 1;
		for (int n2 = 0; n2 < itr; n2++)
		{
			if (n2 != n1) ww = ww * (D(n2) - xx);
		}
		yy = yy + q(n1) * q(n1) * ww;

	}
	zz = (aa - xx) * vv - yy;

	return zz;
}
//=============================================================================
double SVM::NewEnergy(vector <Vector3d>& Basis, MatrixXd C, VectorXd D, double E, double EE)
{
	int itr = Basis.size() - 1;
	
	VectorXd c = VectorXd::Zero(itr + 1);
	VectorXd Overlap = VectorXd::Zero(itr);
	VectorXd q = VectorXd::Zero(itr);
	double aa = 0;
	double NN = 0;

	//solving the equation for the new eigenvalue===============================================  
	int count = 0;
	double e1 = E;
	double e2 = E - abs(0.5 * (E - EE));
	double e3 = E;

	vector<double> overlap_vector(itr);
	for (int k1 = 0; k1 < itr; k1++)
	{
		overlap_vector[k1] = me.overlap(Basis[itr], Basis[k1][0], Basis[k1][1], Basis[k1][2]);
	}
	vector<double> energy_vector(itr);
	for (int k1 = 0; k1 < itr; k1++)
	{
		energy_vector[k1] = me.energy(Basis[itr], Basis[k1][0], Basis[k1][1], Basis[k1][2]);
	}
	//
	for (int k1 = 0; k1 < itr; k1++)
	{
		for (int k2 = 0; k2 < itr; k2++)
		{
			Overlap(k1) = Overlap(k1) + C(k2, k1) * overlap_vector[k2];
		}
	}
	
	c(itr) = 1;
	for (int k1 = 0; k1 < itr; k1++)
	{
		for (int k2 = 0; k2 < itr; k2++)
		{
			c(k1) = c(k1) - Overlap(k2) * C(k1, k2);
		}
	}
	for (int k1 = 0; k1 < itr; k1++)
	{
		for (int k2 = 0; k2 < itr; k2++)
		{

			NN = NN + c(k1) * c(k2) * NM(k1, k2);
		}
	}
	for (int k1 = 0; k1 < itr; k1++)
	{
		NN = NN + 2 * c(k1) * c(itr) * overlap_vector[k1];
	}
	NN = NN + me.overlap(Basis[itr], Basis[itr][0], Basis[itr][1], Basis[itr][2]);
	for (int k1 = 0; k1 < itr + 1; k1++)
	{
		c(k1) = c(k1) / sqrt(NN);
	}
	//------------------------------------------------------
	
	VectorXd HMc=VectorXd::Zero(itr);
	for (int k2 = 0; k2 < itr; k2++) {
		for (int k3 = 0; k3 < itr; k3++) {
			HMc(k2) = HMc(k2)+HM(k2, k3) * c(k3);
		}
	}

	for (int k1 = 0; k1 < itr; k1++)
	{
		for (int k2 = 0; k2 < itr; k2++)
		{
				q(k1) = q(k1) + HMc(k2) * C(k2, k1);
		}
	}
	
	for (int k1 = 0; k1 < itr; k1++)
	{
		for (int k2 = 0; k2 < itr; k2++)
		{
			q(k1) = q(k1) + energy_vector[k2] * c(itr) * C(k2, k1);
		}
	}
	for (int k1 = 0; k1 < itr; k1++)
	{
		for (int k2 = 0; k2 < itr; k2++)
		{
			aa = aa + c(k1) * c(k2) * HM(k1, k2);
		}
	}

	for (int k1 = 0; k1 < itr; k1++)
	{
		aa = aa + 2 * c(k1) * c(itr) * energy_vector[k1];
	}
	aa = aa + c(itr) * c(itr) * me.energy(Basis[itr], Basis[itr][0], Basis[itr][1], Basis[itr][2]);
	
	double Ee1 = EigenValuesEquation(itr, D, q, aa, e1);
	while (count < 101)
	{

		if (Ee1 * EigenValuesEquation(itr, D, q, aa, e2) < 0)  break;
		else
		{
			e1 = e2;
			e2 = e2 - abs(0.5 * (E - EE));
			count++;
		}

	}
	if (count <= 100)
	{
		count = 0;
		while (abs((e1 - e2) / e2) > abs(1e-5 * (E - EE) / EE))
		{
			e3 = (e1 + e2) / 2;
			if (std::signbit(EigenValuesEquation(itr, D, q, aa, e3)) != std::signbit(EigenValuesEquation(itr, D, q, aa, e2)))  e1 = e3;
			else e2 = e3;
			count++;
			if (count > 100)  break;
		}

	}

	return e3;
	//==============================================================================
}


//=============================================================================
MatrixXd SVM::NormMatrix(vector<Vector3d>& Basis)
{
	int itr = Basis.size();
	MatrixXd Norm = MatrixXd::Zero(itr, itr);
	for (int n1 = 0; n1 < itr; n1++)
	{
		for (int n2 = n1; n2 < itr; n2++)
		{

			Norm(n1, n2) = me.overlap(Basis[n1], Basis[n2][0], Basis[n2][1], Basis[n2][2]);
			Norm(n2, n1) = Norm(n1, n2);
		}
	}
	return Norm;
}
//=============================================================================
MatrixXd SVM::HamiltonianMatrix(vector<Vector3d>& Basis)
{
	int itr = Basis.size();
	MatrixXd H = MatrixXd::Zero(itr, itr);
	for (int n1 = 0; n1 < itr; n1++)
	{
		for (int n2 = n1; n2 < itr; n2++)
		{
			H(n1, n2) = me.energy(Basis[n1], Basis[n2][0], Basis[n2][1], Basis[n2][2]);
			H(n2, n1) = H(n1, n2);
		}
	}
	return H;
}
//=============================================================================
void SVM::UpdateNorm(vector<Vector3d>& Basis)
{
	int ncur = Basis.size() - 1;
	for (int n1 = 0; n1 <= ncur; n1++)
	{
		//      cout << Basis.size() << " " << n1 << "\n";
		Nmatrix(n1, ncur) = me.overlap(Basis[n1], Basis[ncur][0], Basis[ncur][1], Basis[ncur][2]);
		Nmatrix(ncur, n1) = Nmatrix(n1, ncur);
	}
}
//=============================================================================
void SVM::UpdateHamiltonian(vector<Vector3d>& Basis)
{
	int ncur = Basis.size() - 1;
	for (int n1 = 0; n1 <= ncur; n1++)
	{
		Hmatrix(n1, ncur) = me.energy(Basis[n1], Basis[ncur][0], Basis[ncur][1], Basis[ncur][2]);
		Hmatrix(ncur, n1) = Hmatrix(n1, ncur);
	}
}