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
int main()
{
    int maxbasis = 100;
    SVM svm;
    Solver solver;
    Vector3d NewState;
    GeneralizedSelfAdjointEigenSolver<MatrixXd> ges;
    vector<Vector3d> Basis;
    
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
    MatrixXd C;
    VectorXd D;
    double E;
    double EE;

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
        clock_t end1 = clock();
        std::cout << "duration = " << (double)(end1 - start1) / CLOCKS_PER_SEC << "sec.\n";

        //Norm matrix
        NM = svm.Nmatrix;
        //Hamiltonian matrix
        HM = svm.Hmatrix;
        solver.Initialize(Basis, C, D, E, EE);
        NewState = solver.GAsolve(Basis, C, D, E, EE);
        

        Basis.push_back(NewState);
        // now the basis size is itr+1

        cout << "itr\t" << itr << "newstate is\t" <<endl<< NewState<<endl;
        cout << "Energy is " << svm.NewEnergy(Basis, C, D, E, EE) << endl;
        EE = E;
        svm.UpdateNorm(Basis);
        svm.UpdateHamiltonian(Basis);
        itr = itr + 1;

    }
    return 0;
}