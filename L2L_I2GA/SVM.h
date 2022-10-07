#pragma once
#define SVM_H
#include "MatrixElement.h"
#include "Eigen/Core"
#include <vector>
using namespace Eigen;
using namespace std;
extern MatrixXd NM;
extern MatrixXd HM;
extern MatrixXd C;
extern VectorXd D;
extern double E;
extern double EE;
extern vector<Vector3d>Basis;
class SVM
{

private:
	MatrixElement me;
	int  kk0;
	double bmin, bmax;

public:
	SVM();

	int CheckOverlap(vector<Vector3d>Basis);
	double NewEnergy();
	MatrixXd NormMatrix(vector<Vector3d>Basis);
	MatrixXd HamiltonianMatrix(vector<Vector3d>Basis);
	MatrixXd Hmatrix, Nmatrix;
	MatrixXd H, Norm;
	void UpdateNorm(vector<Vector3d>Basis);
	void UpdateHamiltonian(vector<Vector3d>Basis);
};
