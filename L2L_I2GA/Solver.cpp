#include "Solver.h"
#include "SpecimenWorkshop.h"
#include <iostream>
#include <vector>

// Current Population
Specimen** Solver::Current_Population;

// Current Selection
Specimen** Solver::Current_Selection;

// Current Proximity
double Solver::Current_Proximity;

// Current Iteration
int Solver::Current_Iteration;

// Initialize the algorithm
void Solver::first_state_Initialize()
{
	SpecimenWorkshop sw;

	// Generate initial population
	Current_Population = sw.first_state_GeneratePopulation();

	// Set Current_Selection to zero for correct work delete[] operator
	Current_Selection = 0;

	Current_Iteration = 0;
}
Vector3d Solver::first_state_GAsolve()
{
	int count = 0;
	// Set Current_Proximity to the biggest value
	Current_Proximity = 1.7976931348623157E+308;
	SpecimenWorkshop sw;
	// Loop while Current_Proximity is not less than the Epsilon
	while (1)
	{
		double Previous_Proximity = Current_Proximity;
		// Delete old selection
		if (Current_Selection != 0)
		{
			for (int i = 0; i < sw.SelectionSize; i++)
			{
				// Calculate Energy for new specimens
				delete Current_Selection[i];
			}

			delete[] Current_Selection;
		}

		// Select the best specimens
		Current_Selection = sw.Select(Current_Population);

		// Calculate proximity for the top-selected (the best) specimen
		Current_Proximity = Current_Selection[0]->Energy;
		Current_Selection[0]->Genes[0];
		Current_Selection[0]->Genes[1];
		Current_Selection[0]->Genes[2];
		// End the calculations if Current_Proximity is less than the Epsilon


		double diff = Previous_Proximity - Current_Proximity;
		if (diff >= 0 && diff < sw.Epsilon)
		{
			count++;
		}


		// Delete old population
		for (int i = 0; i < sw.PopulationSize; i++)
		{
			// Calculate Energy for new specimens
			delete Current_Population[i];
		}

		delete[] Current_Population;
		if (count == 5) {
			Vector3d b;
			double f1 = Current_Selection[0]->Genes[0] * Current_Selection[0]->Genes[0];
			double f2 = Current_Selection[0]->Genes[1] * Current_Selection[0]->Genes[1];
			double f3 = Current_Selection[0]->Genes[2] * Current_Selection[0]->Genes[2];
			b[0] = f1 + 0.25 * f2 + 0.25 * f3;
			b[1] = 0.5 * f2 - 0.5 * f3;
			b[2] = f2 + f3;
			return b;
			break;
		}
		// Generate new population by reproducing specimens from the selection
		Current_Population = sw.first_state_GeneratePopulation(Current_Selection);

		Current_Iteration++;

	}
}
void Solver::Initialize(vector<Vector3d>& Basis, MatrixXd C, VectorXd D, double E, double EE)
{
	SpecimenWorkshop sw;

	// Generate initial population
	Current_Population = sw.GeneratePopulation(Basis, C, D, E, EE);

	// Set Current_Selection to zero for correct work delete[] operator
	Current_Selection = 0;

	Current_Iteration = 0;
}

// Run the algorithm
Vector3d Solver::GAsolve(vector<Vector3d>& Basis, MatrixXd C, VectorXd D, double E, double EE)
{
	int Bsize = Basis.size();
	int count = 0;
	// Set Current_Proximity to the biggest value
	Current_Proximity = 1.7976931348623157E+308;
	SpecimenWorkshop sw;
	// Loop while Current_Proximity is not less than the Epsilon
	while (1)
	{
		double Previous_Proximity = Current_Proximity;
		// Delete old selection
		if (Current_Selection != 0)
		{
			//cout << "ok" << endl;
			for (int i = 0; i < sw.SelectionSize; i++)
			{
				// Calculate Energy for new specimens
				delete Current_Selection[i];
			}

			delete[] Current_Selection;
		}

		// Select the best specimens
		Current_Selection = sw.Select(Current_Population);

		// Calculate proximity for the top-selected (the best) specimen
		Current_Proximity = Current_Selection[0]->Energy;
		Current_Selection[0]->Genes[0];
		Current_Selection[0]->Genes[1];
		Current_Selection[0]->Genes[2];
		// End the calculations if Current_Proximity is less than the Epsilon
		
		
		double diff = Previous_Proximity - Current_Proximity;
		if (diff>=0 && diff< sw.Epsilon)
		{
			count++;
		}
		
		
		// Delete old population
		for (int i = 0; i < sw.PopulationSize; i++)
		{
			// Calculate Energy for new specimens
			delete Current_Population[i];
		}

		delete[] Current_Population;
		if (count == 5) {
			Vector3d b;
			double f1 = Current_Selection[0]->Genes[0] * Current_Selection[0]->Genes[0];
			double f2 = Current_Selection[0]->Genes[1] * Current_Selection[0]->Genes[1];
			double f3 = Current_Selection[0]->Genes[2] * Current_Selection[0]->Genes[2];
			b[0] = f1 + 0.25 * f2 + 0.25 * f3;
			b[1] = 0.5 * f2 - 0.5 * f3;
			b[2] = f2 + f3;
			return b;
			break;
		}
		// Generate new population by reproducing specimens from the selection
		Current_Population = sw.GeneratePopulation(Current_Selection, Basis, C, D, E, EE);

		Current_Iteration++;
		
	}
}
