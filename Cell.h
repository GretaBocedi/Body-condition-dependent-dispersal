/*------------------------------------------------------------------------------
Body condition-dependent dispersal: Cell

Implements the Cell class
Each cell in the landscape contains an instance of this class.

The function outPop is defined in BCDispersal.cpp

THE JOINT EVOLUTION OF DENSITY- AND BODY CONDITION-DEPENDENT DISPERSAL
Baines, C.B., Travis, J.M.J, McCauley, S.J. & Bocedi, G. - submitted

Author: Greta Bocedi

Last updated: 29/01/2020 by G. Bocedi
--------------------------------------------------------------------------------*/

#pragma once

#include <vector>
#include "Individual.h"

class Cell
{
public:
	Cell(int, int);
	~Cell();
	int x, y; // cell coordinates
	int N, Nfem, Nmale; //number of individual in the population
	int Noffs, Foffs, Moffs; //number of offspring in the population
	int emigrants; //nr. of individuals emigrated from the cell
	double Res; //amount of resources
	double prop_emigrants;
	double eps; //environmental noise

	double sum_BC; //sum of all alive individuals body condition
	
	std::vector<Individual> females, males, Jfemales, Jmales, survJfem, survJmales;

	void outPop(void); //function to output population's attributes
	void DeleteAdults(void);
};

