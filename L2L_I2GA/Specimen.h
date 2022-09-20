#ifndef H_SPECIMEN
#define H_SPECIMEN
#include "MatrixElement.h"
#include <vector>
#include "Eigen/Core"
#include "SVM.h"

using namespace Eigen;
using namespace std;
// Class that represents the Specimen
class Specimen
{
private:
public:
	// The genes
	double* Genes;

	// Energy to the solution
	double Energy;

	// Constructor that creates the Specimen
	// instance with initial genes
	Specimen();

	// Destructor. Frees memory for the Genes
	~Specimen();
	// Clone the Specimen for simple memory management
	Specimen* Clone();

	// Calculates Energy of this instance to the solution
	// Contains the equation formula
	void CalculateEnergy(vector<Vector3d>& Basis, MatrixXd C, VectorXd D, double E, double EE);
	void first_state_CalculateEnergy();
};

#endif
