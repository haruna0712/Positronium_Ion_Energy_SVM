#ifndef H_SOLVER
#define H_SOLVER

#include "Specimen.h"
#include <vector>
// Class that implements the Genetic Algorithm
// at the heighest absraction level
class Solver
{
public:
	// Current Population
	static Specimen** Current_Population;

	// Current Selection
	static Specimen** Current_Selection;

	// Current Proximity
	static double Current_Proximity;

	// Current Iteration
	static int Current_Iteration;

	void first_state_Initialize();

	static void Initialize(vector<Vector3d>Basis);
	
	Vector3d GAsolve(vector<Vector3d>Basis);

	Vector3d first_state_GAsolve();
};

#endif
