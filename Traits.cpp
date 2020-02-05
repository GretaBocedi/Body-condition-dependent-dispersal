#include "Traits.h"

//Calss constructor
Traits::Traits()
{
	alpha = NULL;
	beta = NULL;
	D0 = NULL;
	alphaD = NULL;
	betaD = NULL;
	gamma = NULL;

	g_alpha = NULL;
	g_beta = NULL;
	g_D0 = NULL;
	g_alphaD = NULL;
	g_betaD = NULL; 
	g_gamma = NULL;

	p_alpha = NULL;
	p_beta = NULL;
	p_D0 = NULL;
	p_alphaD = NULL;
	p_betaD = NULL;
	p_gamma = NULL;
}

//--------------------------------------------------------
Traits::~Traits()
{
}
//--------------------------------------------------------
void Traits::deleteTraits(void){
	if (alpha != NULL){ delete[] alpha; alpha = NULL; }
	if (beta != NULL){ delete[] beta; beta = NULL; }
	if (D0 != NULL){ delete[] D0; D0 = NULL; }
	if (alphaD != NULL){ delete[] alphaD; alphaD = NULL; }
	if (betaD != NULL){ delete[] betaD; betaD = NULL; }
	if (gamma != NULL){ delete[] gamma; gamma = NULL; }

	if (g_alpha != NULL){ delete g_alpha; g_alpha = NULL; }
	if (g_beta != NULL){ delete g_beta; g_beta = NULL; }
	if (g_D0 != NULL){ delete g_D0; g_D0 = NULL; }
	if (g_alphaD != NULL){ delete g_alphaD; g_alphaD = NULL; }
	if (g_betaD != NULL){ delete g_betaD; g_betaD = NULL; }
	if (g_gamma != NULL){ delete g_gamma; g_gamma = NULL; }

	if (p_alpha != NULL){ delete p_alpha; g_alpha = NULL; }
	if (p_beta != NULL){ delete p_beta; p_beta = NULL; }
	if (p_D0 != NULL){ delete p_D0; p_D0 = NULL; }
	if (p_alphaD != NULL){ delete p_alphaD; p_alphaD = NULL; }
	if (p_betaD != NULL){ delete p_betaD; p_betaD = NULL; }
	if (p_gamma != NULL){ delete p_gamma; p_gamma = NULL; }
}
