/*------------------------------------------------------------------------------
Body condition-dependent dispersal: Cell

Implements the Traits class
Each individual contains an instance of this class.

THE JOINT EVOLUTION OF DENSITY- AND BODY CONDITION-DEPENDENT DISPERSAL
Baines, C.B., Travis, J.M.J, McCauley, S.J. & Bocedi, G. - submitted

Author: Greta Bocedi

Last updated: 29/01/2020 by G. Bocedi
--------------------------------------------------------------------------------*/

#pragma once
#include <stdio.h>


class Traits
{
public:
	Traits();
	~Traits();
	//alleles
	double *alpha; //slope of the dispersal reaction norm to body-condition
	double *beta;  //intercept of the dispersal reaction norm to body-condition
	double *D0; //maximum dispersal probability
	double *alphaD; //slope of the dispersal reaction norm to density
	double *betaD; //intercept of the dispersal reaction norm to density
	double *gamma; //interaction term between body condition and density
	//genotypes
	double *g_alpha, *g_beta, *g_D0, *g_alphaD, *g_betaD, *g_gamma;
	//phenotypes
	double *p_alpha, *p_beta, *p_D0, *p_alphaD, *p_betaD, *p_gamma;

	void deleteTraits(void);
};

