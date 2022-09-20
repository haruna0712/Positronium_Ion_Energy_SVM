#include "Specimen.h"
#include <math.h>
#include <cmath>
#include <time.h>
#include "SVM.h"
#include <vector>
#include "MatrixElement.h"
#include<iostream>

MatrixElement me;
SVM svm;
Specimen::Specimen()
{
	// Assume one double value as the genes
	Genes = new double[3];

	// Set up initial value to the genes
	//Genes[0] is the coefficient of (r1-r2)^2,Genes[1] is the coefficient of (r2-r3)^2,Genes[2] is the coefficient of (r3-r1)^2,
	Genes[0] = 10;
	Genes[1] = 10;
	Genes[2] = 10;
}

// Destructor. Frees memory for the Genes
Specimen::~Specimen()
{
	delete[] Genes;
}

// Clone the Specimen for simple memory management
Specimen* Specimen::Clone()
{
	Specimen* s = new Specimen();
	s->Genes[0] = Genes[0];
	s->Genes[1] = Genes[1];
	s->Genes[2] = Genes[2];
	s->Energy = Energy;

	return s;
}
void Specimen::CalculateEnergy(vector<Vector3d>& Basis, MatrixXd C, VectorXd D, double E, double EE)
{
	int Bsize = Basis.size();
	Vector3d NewState;
	double f1 = Genes[0] * Genes[0];
	double f2 = Genes[1] * Genes[1];
	double f3 = Genes[2] * Genes[2];

	//coodinate transformation from descartes coordinate (r1,r1,r3) to relative (Jacobi) coordinate (r12,rho)
	NewState[0] = f1 + 0.25 * f2 + 0.25 * f3;
	NewState[1] = 0.5 * (f2 - f3);
	NewState[2] = f2 + f3;
	Basis.push_back(NewState);
	
	double NewE;
	NewE = svm.NewEnergy(Basis, C, D, E, EE);

	Energy = NewE;
	Basis.resize(Bsize);
}
void Specimen::first_state_CalculateEnergy()
{
	Vector3d NewState;
	double f1 = Genes[0] * Genes[0];
	double f2 = Genes[1] * Genes[1];
	double f3 = Genes[2] * Genes[2];

	//coodinate transformation from single particle coordinate (r1,r1,r3) to relative (Jacobi) coordinate (r12,rho)

	double b11 = f1 + 0.25 * f2 + 0.25 * f3;
	double b12 = 0.5 * f2 - 0.5 * f3;
	double b22 = f2 + f3;

	NewState[0] = b11;
	NewState[1] = b12;
	NewState[2] = b22;

	Energy = me.energy(NewState, b11,b12,b22) / me.overlap(NewState, b11, b12, b22);
}