#include "SVM.h"
#include "MatrixElement.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include "Eigen/Core"
#include "Eigen/Eigen"
#include <iomanip>
#include <string>
#include <stdlib.h> 
#include <ctime> 
#include <cstdio>
#include <cstdlib>
#include <stdio.h>
#include <time.h>
#include "Solver.h"
using namespace Eigen;
using namespace std;

MatrixXd NM;
MatrixXd HM;
MatrixXd C;
VectorXd D;
double E;
double EE;
vector<Vector3d>Basis;
int main()
{
    int maxbasis = 200;
    SVM svm;
    Solver solver;
    Vector3d NewState;
    GeneralizedSelfAdjointEigenSolver<MatrixXd> ges;
    
    solver.first_state_Initialize();
    NewState = solver.first_state_GAsolve();
    cout << "first new state is  " << endl;
    cout << NewState << "\n";

    /* start SVM iterations */
    Basis.push_back(NewState);
    svm.UpdateNorm(Basis);
    svm.UpdateHamiltonian(Basis);

    MatrixXd Norm;
    MatrixXd H;

    cout << "\t Start SVM iters\n\n";
    int itr = 1;
    clock_t start1 = clock();
    //-------------------------------------------------
    while (itr < maxbasis)
    {

        Norm = svm.NormMatrix(Basis);
        H = svm.HamiltonianMatrix(Basis);
        ges.compute(H, Norm);
        C = ges.eigenvectors();
        D = ges.eigenvalues();
        E = D.minCoeff();
        if (itr == 1)  EE = E + abs(E / 2);

        printf("\t iter = %4d     E = %14.8f  \n", itr, E);
        //Norm matrix
        NM = svm.Nmatrix;
        //Hamiltonian matrix
        HM = svm.Hmatrix;
        solver.Initialize(Basis);
        NewState = solver.GAsolve(Basis);

        Basis.push_back(NewState);
        cout << "Energy is " << svm.NewEnergy() << endl;
        EE = E;
        svm.UpdateNorm(Basis);
        svm.UpdateHamiltonian(Basis);
        itr = itr + 1;

    }
    clock_t end1 = clock();
    std::cout << "update duration = " << (double)(end1 - start1) / CLOCKS_PER_SEC << "sec.\n";
    return 0;
}