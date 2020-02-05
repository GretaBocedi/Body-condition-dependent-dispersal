/*------------------------------------------------------------------------------
Body condition-dependent dispersal: Parameters

Implements the Parameters class

The function outPop is defined in BCDispersal.cpp

THE JOINT EVOLUTION OF DENSITY- AND BODY CONDITION-DEPENDENT DISPERSAL
Baines, C.B., Travis, J.M.J, McCauley, S.J. & Bocedi, G. - submitted

Author: Greta Bocedi

Last updated: 29/01/2020 by G. Bocedi
--------------------------------------------------------------------------------*/
#pragma once
class Parameters
{
public:
	Parameters();
	~Parameters();
	//Simulation
	int SimNr; //Simulation number
	int rep; //replicates
	int gen; //generations
	int ind_int; //generations interval for individual outputs
	int pop_int; //generations interval for population outputs
	//Landscape
	int x_max;
	int y_max;
	int max_seedX; //max x for initialisation
	int max_seedY; //max y for initialisation
	double minR; //minimum cell's amount of resources
	double maxR; //maximum cell's amount of resources
	double minpar; //minimum resource parcel size
	double maxpar; //maximum resource parcel size

	//Traits-----------------------------------------------
	//Parameters for BC-dependent dispersal
	double alpha_Imean; // initial genotypic mean for alpha
	double beta_Imean; // initial genotypic mean for beta
	double D0_Imean; // initial genotypic mean for D0
	double alpha_Isd; // initial genotypic standard deviation for alpha
	double beta_Isd; // initial genotypic standard deviation for beta
	double D0_Isd; // initial genotypic standard deviation for D0
	//Parameters for density-dependent dispersal (use D0 as max dispersal distance)
	double alphaD_Imean; // initial genotypic mean for alphaD
	double betaD_Imean; // initial genotypic mean for betaD
	double alphaD_Isd; // initial genotypic standard deviation for alphaD
	double betaD_Isd; // initial genotypic standard deviation for betaD
	//Interaction between BC and density 
	double gamma_Imean; // initial genotypic mean for gamma
	double gamma_Isd; // initial genotypic standard deviation for gamma
	//Density and BC independent dispersal: use D0.
	//Mutations-----------------------------------------------
	double mu; //haploid per allele mutation rate for each trait
	//Parameters for BC-dependent dispersal
	double mu_mean_alpha; //mean of the distribution of mutational effects for alpha
	double mu_sd_alpha; //standard deviation of the distribution of mutational effects for alpha
	double mu_mean_beta; //mean of the distribution of mutational effects for beta
	double mu_sd_beta; //standard deviation of the distribution of mutational effects for beta
	double mu_mean_D0; //mean of the distribution of mutational effects for D0
	double mu_sd_D0; //standard deviation of the distribution of mutational effects for D0
	//Parameters for density-dependent dispersal
	double mu_mean_alphaD; //mean of the distribution of mutational effects for alphaD
	double mu_sd_alphaD; //standard deviation of the distribution of mutational effects for alphaD
	double mu_mean_betaD; //mean of the distribution of mutational effects for betaD
	double mu_sd_betaD; //standard deviation of the distribution of mutational effects for betaD
	//Interaction between BC and density 
	double mu_mean_gamma; //mean of the distribution of mutational effects for gamma
	double mu_sd_gamma; //standard deviation of the distribution of mutational effects for gamma

	//Resource-dependent body condition
	double alpha_res; //slope of logistic function inflection point
	double beta_res; //inflection point
	//Reproduction
	bool BCfecundity; //BC-dependent fecundity
	double max_fec; //maximum fecundity for body-condition = 1
	//Dispersal
	int disp_function;	// 0 = linear; 
						// 1 = body condition dep. (3-parameters logistic); 
						// 2 = body-condition and density dependent disp. (5 parameters)
						// 3 = density-dependent (3-parameters logistic); 
						// 4 = condition- (body or density) -independent dispersal
	int disp_type; // 0 = nearest_neighbour; 1 = global dispersal; 
	double mean_dist; //mean dispersal distance in the case of dispersal kernel (in cell unit)
	//Survival dispersal cost
	int disp_mort_type; // 0 = constant (disp_mort); // 1 = BC-dependent (beta_surv, alpha_surv)
	double disp_mort; //constant dispersal mortality
	double beta_surv; //inflection point BC-dependent dispersal survival
	double alpha_surv; //slope BC-dependent dispersal survival
	//Energetic dispersal cost
	double disp_cost; //dispersal cost paid as reduction in BC (hence fecundity for females)
	//Apply survival after dispersal
	bool survival; //yes/no
	//Environmental stochasticity
	double ac; //temporal autocorrelation (-9 = no env. stochasticity)
	double std; //standard deviation of env. fluctuations
	double local_ext; //local extinction probability

};

