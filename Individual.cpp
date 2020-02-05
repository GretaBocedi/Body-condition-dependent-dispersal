#include "Individual.h"

// Random numbers generators and seed 
std::random_device rd2;
std::mt19937 gen(rd2());

// Class constructor
Individual::Individual(bool sx, int xx, int yy)
{
	sex = sx;
	alive = true;
	dispersed = false;
	x = xx;
	y = yy;
	new_x = -9;
	new_y = -9;
	n_par = 0;
	resources = 0.0;
	BC = 0.0;
	d = 0.0;
}
//-----------------------------------------------------------
Individual::~Individual()
{
}
//---------------------------------------------------------------------------
//Initialise individual dispersal traits (alleles, genotypes and phenotypes)
void Individual::initialise(int type, std::normal_distribution<> distrAlpha, std::normal_distribution<> distrBeta, std::normal_distribution<> distrD0,
	std::normal_distribution<> distrAlphaD, std::normal_distribution<> distrBetaD, std::normal_distribution<> distrGamma) {
	
	switch (type) {
	case 0: //linear function
			//alpha
		traits.alpha = new double[2];
		traits.g_alpha = new double(0.0);
		traits.p_alpha = new double(0.0);
		//beta
		traits.beta = new double[2];
		traits.g_beta = new double(0.0);
		traits.p_beta = new double(0.0);
		//alleles
		traits.alpha[0] = distrAlpha(gen); traits.alpha[1] = distrAlpha(gen);
		traits.beta[0] = distrBeta(gen); traits.beta[1] = distrBeta(gen);
		//Genotypes (sum of alleles)
		*traits.g_alpha = traits.alpha[0] + traits.alpha[1];
		*traits.g_beta = traits.beta[0] + traits.beta[1];
		//Phenotypes
		*traits.p_alpha = *traits.g_alpha;
		*traits.p_beta = *traits.g_beta;
		break;
	case 1: //body condition dep. (3-parameters logistic); 
			//alpha
		traits.alpha = new double[2];
		traits.g_alpha = new double(0.0);
		traits.p_alpha = new double(0.0);
		//beta
		traits.beta = new double[2];
		traits.g_beta = new double(0.0);
		traits.p_beta = new double(0.0);
		//D0
		traits.D0 = new double[2];
		traits.g_D0 = new double(0.0);
		traits.p_D0 = new double(0.0);
		//alleles
		traits.alpha[0] = distrAlpha(gen); traits.alpha[1] = distrAlpha(gen);
		traits.beta[0] = distrBeta(gen); traits.beta[1] = distrBeta(gen);
		traits.D0[0] = distrD0(gen); traits.D0[1] = distrD0(gen);
		//Genotypes (sum of alleles)
		*traits.g_alpha = traits.alpha[0] + traits.alpha[1];
		*traits.g_beta = traits.beta[0] + traits.beta[1];
		*traits.g_D0 = traits.D0[0] + traits.D0[1];
		//Phenotypes
		*traits.p_alpha = *traits.g_alpha;
		*traits.p_beta = *traits.g_beta;
		*traits.p_D0 = *traits.g_D0;
		if (*traits.p_D0 < 0.0) *traits.p_D0 = 0.0; // Constrain D0 between 0 and 1
		if (*traits.p_D0 > 1.0) *traits.p_D0 = 1.0;
		break;
	case 2: //body - condition and density dependent disp. (5 parameters)
			//alpha
		traits.alpha = new double[2];
		traits.g_alpha = new double(0.0);
		traits.p_alpha = new double(0.0);
		//beta
		traits.beta = new double[2];
		traits.g_beta = new double(0.0);
		traits.p_beta = new double(0.0);
		//alphaD
		traits.alphaD = new double[2];
		traits.g_alphaD = new double(0.0);
		traits.p_alphaD = new double(0.0);
		//betaD
		traits.betaD = new double[2];
		traits.g_betaD = new double(0.0);
		traits.p_betaD = new double(0.0);
		//interaction
		traits.gamma = new double[2];
		traits.g_gamma = new double(0.0);
		traits.p_gamma = new double(0.0);
		//alleles
		traits.alpha[0] = distrAlpha(gen); traits.alpha[1] = distrAlpha(gen);
		traits.beta[0] = distrBeta(gen); traits.beta[1] = distrBeta(gen);
		traits.alphaD[0] = distrAlphaD(gen); traits.alphaD[1] = distrAlphaD(gen);
		traits.betaD[0] = distrBetaD(gen); traits.betaD[1] = distrBetaD(gen);
		traits.gamma[0] = distrGamma(gen); traits.gamma[1] = distrGamma(gen);
		//Genotypes (sum of alleles)
		*traits.g_alpha = traits.alpha[0] + traits.alpha[1];
		*traits.g_beta = traits.beta[0] + traits.beta[1];
		*traits.g_alphaD = traits.alphaD[0] + traits.alphaD[1];
		*traits.g_betaD = traits.betaD[0] + traits.betaD[1];
		*traits.g_gamma = traits.gamma[0] + traits.gamma[1];
		//Phenotypes
		*traits.p_alpha = *traits.g_alpha;
		*traits.p_beta = *traits.g_beta;
		*traits.p_alphaD = *traits.g_alphaD;
		*traits.p_betaD = *traits.g_betaD;
		*traits.p_gamma = *traits.g_gamma;
		break;
	case 3: //density-dependent (3-parameters logistic); 
			//alphaD
		traits.alphaD = new double[2];
		traits.g_alphaD = new double(0.0);
		traits.p_alphaD = new double(0.0);
		//betaD
		traits.betaD = new double[2];
		traits.g_betaD = new double(0.0);
		traits.p_betaD = new double(0.0);
		//D0
		traits.D0 = new double[2];
		traits.g_D0 = new double(0.0);
		traits.p_D0 = new double(0.0);
		//alleles
		traits.alphaD[0] = distrAlphaD(gen); traits.alphaD[1] = distrAlphaD(gen);
		traits.betaD[0] = distrBetaD(gen); traits.betaD[1] = distrBetaD(gen);
		traits.D0[0] = distrD0(gen); traits.D0[1] = distrD0(gen);
		//Genotypes (sum of alleles)
		*traits.g_alphaD = traits.alphaD[0] + traits.alphaD[1];
		*traits.g_betaD = traits.betaD[0] + traits.betaD[1];
		*traits.g_D0 = traits.D0[0] + traits.D0[1];
		//Phenotypes
		*traits.p_alphaD = *traits.g_alphaD;
		*traits.p_betaD = *traits.g_betaD;
		*traits.p_D0 = *traits.g_D0;
		if (*traits.p_D0 < 0.0) *traits.p_D0 = 0.0; // Constrain D0 between 0 and 1
		if (*traits.p_D0 > 1.0) *traits.p_D0 = 1.0;
		break;
	case 4: //condition independent dispersal
			//D0
		traits.D0 = new double[2];
		traits.g_D0 = new double(0.0);
		traits.p_D0 = new double(0.0);
		//alleles
		traits.D0[0] = distrD0(gen); traits.D0[1] = distrD0(gen);
		//Genotypes (sum of alleles)
		*traits.g_D0 = traits.D0[0] + traits.D0[1];
		//Phenotypes
		*traits.p_D0 = *traits.g_D0;
		if (*traits.p_D0 < 0.0) *traits.p_D0 = 0.0; // Constrain D0 between 0 and 1
		if (*traits.p_D0 > 1.0) *traits.p_D0 = 1.0;
		break;
	}

}
//-----------------------------------------------------------
//determines the individual body condition
void Individual::bodyCond(double a, double b){
	if (resources == 0.0) alive = false;
	else BC = 1.0 / (1.0 + exp(-(resources - b) * a)); 
}
//-----------------------------------------------------------
//Calculates the individual's dispersal probability
void Individual::dispersalP(int funct, double den){
	switch (funct){
	case 0:
		d = (*traits.p_beta) + (*traits.p_alpha) * BC; //From Bonte & De La Pena, 2009 - J.Evol.Biol.
		break;
	case 1: //BC-dependent dispersal: 3-parameters logistic function
		d = (*traits.p_D0) / (1.0 + exp(-(BC - (*traits.p_beta))*(*traits.p_alpha)));
		break;
	case 2: //body condition and density-dependent dispersal
		d = 1.0 / (1.0 + exp(-(den - (*traits.p_betaD))*(*traits.p_alphaD) - (BC - (*traits.p_beta))*(*traits.p_alpha) - (*traits.p_gamma)*den*BC));
		break;
	case 3: //density-dependent dispersal: 3-parameters logistic function
		d = (*traits.p_D0) / (1.0 + exp(-(den - (*traits.p_betaD))*(*traits.p_alphaD)));
		break;
	case 4: //condition independent dispersal
		d = *traits.p_D0;
		break;
	}
}
//-----------------------------------------------------------
void Individual::deleteInd(void){
	traits.deleteTraits();
}