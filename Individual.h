/*------------------------------------------------------------------------------
Body condition-dependent dispersal: Individual

Implements the Individual class

The function disperse, mutation and outInd are defined in BCDispersal.cpp

THE JOINT EVOLUTION OF DENSITY- AND BODY CONDITION-DEPENDENT DISPERSAL
Baines, C.B., Travis, J.M.J, McCauley, S.J. & Bocedi, G. - submitted

Author: Greta Bocedi

Last updated: 29/01/2020 by G. Bocedi
--------------------------------------------------------------------------------*/


#pragma once
#include <math.h>
#include <random>
#include "Traits.h"


class Individual
{
public:
	Individual(bool, int, int);
	~Individual();

	bool sex;
	bool alive;
	bool dispersed;
	int x, y; //individual's coordinates
	int new_x, new_y; //individual's coordinates after it dispersed
	int n_par; //nr. of resouce parcels acquired
	double resources; //amount of resources available to the individual
	double BC; //individual's body condition
	double d; //dispersal probability

	Traits traits;

	void initialise(int, std::normal_distribution<>, std::normal_distribution<>, std::normal_distribution<>,
		std::normal_distribution<>, std::normal_distribution<>, std::normal_distribution<>); //function to initialize the individual
	void bodyCond(double a, double b); //determines the individual body condition
	void mutation(int trait, int allele);
	void dispersalP(int funct, double den); //Calculates the individual's dispersal probability
	void disperse(int max_x, int max_y, int type);//calculates the individual's new location
	void outInd(void); //individual output
	void deleteInd(void);

};

