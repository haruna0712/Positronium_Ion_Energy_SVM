#include "MatrixElement.h"
#include <cmath>
#include <iostream>
#include "Eigen/Core"
#include <vector>
using namespace Eigen;
using namespace std;

//=============================================================================
MatrixElement::MatrixElement()
{
}
// the ground state of positronium ion is singlet as for two electrons , so the wave function must be symmetrized in real space. 
// overlap_single means overlap to one component of the symmetrized wave function.
// A is a vector containing the matrix element of A (coefficient of the Gaussian bra state)
// b11, b12,b22 are the matrix element of symmetry matrix B, which is the coefficient of Gaussian ket state.
double MatrixElement::overlap_single(Vector3d A, double b11, double b12, double b22)
{
	double a11 = A[0];
	double a12 = A[1];
	double a22 = A[2];
	double overlap;
	
	double detAplusB = (a11 + b11) * (a22 + b22) - (a12 + b12) * (a12 + b12);
	overlap = pow(detAplusB,-1.5);

	return overlap;
}
double MatrixElement::overlap(Vector3d A, double b11, double b12, double b22)
{
	double TBT11 = b11;
	double TBT12 = -b12;
	double TBT22 = b22;

	return  overlap_single(A,b11,b12,b22)+overlap_single(A, TBT11, TBT12, TBT22);
}

//==============================================================================================
double MatrixElement::energy_single(Vector3d A, double b11,double b12,double b22) {

	double a11 = A[0];
	double a12 = A[1];
	double a22 = A[2];
	double detA;
	double detB;
	double detAplusB;

	detA = a11 * a22 - a12 * a12;
	detB = b11 * b22 - b12 * b12;
	double a12b12 = a12 * b12;

	//determinant of A+B
	detAplusB = (a11 + b11) * (a22 + b22) - (a12 + b12) * (a12 + b12);

	// matrix element representation of kinetic energy= (A*(A+B).inverse()*B*lamda).trace()
	double KinCoeff = (detA * (3 * b11 + 2.25 * b22) + detB * (3 * a11 + 2.25 * a22));
	double KinEnergy = KinCoeff /pow(detAplusB,2.5);
	
	c12_inv = a22 + b22;
	c23_inv = (a22 + b22) * 0.25 + (a12 + b12) + (a11+b11);
	c13_inv = (a22 + b22) * 0.25 - (a12 + b12) + (a11+b11);

	double PotCoeff = (sqrt(1.0 / c12_inv) - sqrt(1.0 / c23_inv) - sqrt(1.0 / c13_inv)) * sqrt(2.0 / pi);

	double PotEnergy = PotCoeff/detAplusB ;
	return KinEnergy + PotEnergy;


}
double MatrixElement::energy(Vector3d A, double b11,double b12,double b22) {
	double E1 = energy_single(A, b11,b12,b22);

	double TBT11 = b11;
	double TBT12 = -b12;
	double TBT22 = b22;

	double E12 = energy_single(A, TBT11,TBT12,TBT22);

	return E1+E12;
}