#include "SpecimenWorkshop.h"
#include <vector>
#include <cmath>
#include <fstream>
#include <random>
#include <algorithm>
#include<iostream>

std::random_device rd;
std::default_random_engine eng(rd());
std::uniform_real_distribution<double> distr(0.0, 1.0);

// Number of specimens in the population
int SpecimenWorkshop::PopulationSize=20;

// Number of specimens in the selection 
int SpecimenWorkshop::SelectionSize=7;


// Likelyhood of mutation
double SpecimenWorkshop::MutationLikelyhoodPercent=0.2;

// Maximal Energy that is considered 
// for the solution found
double SpecimenWorkshop::Epsilon=0.01;

// Generate initial population
Specimen** SpecimenWorkshop::GeneratePopulation(vector<Vector3d>Basis)
{
	// Creates array representing the population
	Specimen** p = new Specimen*[PopulationSize];

	for (int i = 0; i < PopulationSize; i++)
	{
		p[i] = new Specimen();
		Mutate(p[i]);

		// Calculate Energy for new specimens
		p[i]->CalculateEnergy();
	}
	return p;
}
Specimen** SpecimenWorkshop::first_state_GeneratePopulation()
{
	// Creates array representing the population
	Specimen** p = new Specimen * [PopulationSize];

	// Creates specimens
	// Mutation of all specimens in initial population
	// increases variance that increases chance to
	// get better instance.

	for (int i = 0; i < PopulationSize; i++)
	{
		p[i] = new Specimen();
		Mutate(p[i]);

		// Calculate Energy for new specimens
		p[i]->first_state_CalculateEnergy();
	}
	return p;
}

// Generate population by reproduction of selection
Specimen** SpecimenWorkshop::GeneratePopulation(Specimen** selection, vector<Vector3d>Basis)
{
	
	// Creates array representing the population
	Specimen** p = new Specimen*[PopulationSize];

	// Copy instances from the selection to keep them
	// in new generation of population
	for (int i = 0; i < SelectionSize; i++)
	{
		p[i] = selection[i]->Clone();
	}

	// Creates new specimens by reproducing two parents
	// Parents are selected randomly from the selection.
	int parent1_index;
	int parent2_index;

	for (int child_index = 0; child_index < PopulationSize; child_index++)
	{
		// Slect two parents randomly in way
		// they are different instances
			parent1_index = rand() % SelectionSize;
			parent2_index = rand() % SelectionSize;

		// Creates new specimen
		p[child_index] = ReproduceNew(selection[parent1_index], selection[parent2_index],Basis);

	
	}
	return p;
}
Specimen** SpecimenWorkshop::first_state_GeneratePopulation(Specimen** selection)
{
	Specimen** p = new Specimen * [PopulationSize];

	// Copy instances from the selection to keep them
	// in new generation of population
	for (int i = 0; i < SelectionSize; i++)
	{
		p[i] = selection[i]->Clone();
	}

	// Creates new specimens by reproducing two parents
	// Parents are selected randomly from the selection.
	int parent1_index;
	int parent2_index;

	for (int child_index = 0; child_index < PopulationSize; child_index++)
	{
		// Slect two parents randomly in way
		// they are different instances
			parent1_index = rand() % SelectionSize;
			parent2_index = rand() % SelectionSize;

		// Creates new specimen
		p[child_index] = first_state_ReproduceNew(selection[parent1_index], selection[parent2_index]);


	}
	return p;
}
// Reproduce new specimen on base of two parents
Specimen* SpecimenWorkshop::ReproduceNew(Specimen* a, Specimen* b, vector<Vector3d>Basis)
{
	Specimen* s = new Specimen();
	// Inherit genes as the average on the parents' genes
	if (distr(eng)<0.5)
	{
		s->Genes[0] = a->Genes[0];
	}
	else s->Genes[0] = b->Genes[0];
	if (distr(eng)<0.5)
	{
		s->Genes[1] = a->Genes[1];
	}
	else s->Genes[1] = b->Genes[1];
	if (distr(eng)<0.5)
	{
		s->Genes[2] = a->Genes[2];
	}
	else s->Genes[2] = b->Genes[2];

	if (distr(eng) <= MutationLikelyhoodPercent)
	{
		Mutate(s);
	}

	// Calculate Energy for new specimen
	s->CalculateEnergy();

	return s;
}
Specimen* SpecimenWorkshop::first_state_ReproduceNew(Specimen* a, Specimen* b)
{
	Specimen* s = new Specimen();
	// Inherit genes as the average on the parents' genes
	if (distr(eng) < 0.5)
	{
		s->Genes[0] = a->Genes[0];
	}
	else s->Genes[0] = b->Genes[0];
	if (distr(eng) < 0.5)
	{
		s->Genes[1] = a->Genes[1];
	}
	else s->Genes[1] = b->Genes[1];
	if (distr(eng) < 0.5)
	{
		s->Genes[2] = a->Genes[2];
	}
	else s->Genes[2] = b->Genes[2];

	if (distr(eng) <= MutationLikelyhoodPercent)
	{
		Mutate(s);
	}

	// Calculate Energy for new specimen
	s->first_state_CalculateEnergy();

	return s;
}
// Select the best specimens from the population
Specimen** SpecimenWorkshop::Select(Specimen** population)
{
	// Sort population by increasing the Energy
	// The best specimens are moving to start of the array
	Sort(population);

	// Create set of selected specimens
	Specimen** selected = new Specimen*[SelectionSize];

	// Copy best specimens into the selection
	for (int i = 0; i < SelectionSize; i++)
	{
		selected[i] = population[i]->Clone();
	}

	return selected;
}

// Sort the population
void SpecimenWorkshop::Sort(Specimen** population)
{
	std::sort(population, population+PopulationSize,
		[](const auto& lhs, const auto& rhs) {
			return lhs->Energy < rhs->Energy;
		}
	);
}

// Mutate the specimen
void SpecimenWorkshop::Mutate(Specimen* sp)
{
	// Calculate new gene
	sp->Genes[0] = sp->Genes[0] * (distr(eng));
	sp->Genes[1] = sp->Genes[1] * (distr(eng));
	sp->Genes[2] = sp->Genes[2] * (distr(eng));
}
