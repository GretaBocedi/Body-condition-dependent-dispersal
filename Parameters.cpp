#include "Parameters.h"

//Simulation parameters
Parameters::Parameters()
{
	//Simulation
	SimNr = 322;
	rep = 20; //replicates
	gen = 100001; //generations
	ind_int = 100000; //generations interval for individual outputs
	pop_int = 5000; //generations interval for population outputs
	//Landscape
	x_max = 10;
	y_max = 10;
	max_seedX = x_max; //max x for initialisation
	max_seedY = y_max; //max y for initialisation
	minR = 0.0; //minimum cell's amount of resources
	maxR = 100.0; //maximum cell's amount of resources
	minpar = 0.05; //minimum resource parcel size
	maxpar = 0.2; //maximum resource parcel size
	//Traits-----------------------------------------------
	//Parameters for BC-dependent dispersal
	alpha_Imean = 0.0; // initial genotypic mean for alpha
	beta_Imean = 0.5; // initial genotypic mean for beta
	D0_Imean = 0.5; // initial genotypic mean for D0
	alpha_Isd = 1.0; // initial genotypic standard deviation for alpha
	beta_Isd = 0.1; // initial genotypic standard deviation for beta
	D0_Isd = 0.1; // initial genotypic standard deviation for D0
	//Parameters for density-dependent dispersal
	alphaD_Imean = 0.0; // initial genotypic mean for alphaD
	betaD_Imean = 0.5; // initial genotypic mean for betaD
	alphaD_Isd = 1.0; // initial genotypic standard deviation for alphaD
	betaD_Isd = 0.1; // initial genotypic standard deviation for betaD
	//Interaction between BC and density 
	gamma_Imean = 0.0; // initial genotypic mean for gamma
	gamma_Isd = 0.1; // initial genotypic standard deviation for gamma
	//Density and BC independent dispersal: use D0.
	//Mutations-----------------------------------------------
	mu = 0.01; //haploid per allele mutation rate for each trait
	//Parameters for BC-dependent dispersal
	mu_mean_alpha = 0.0; //mean of the distribution of mutational effects for alpha
	mu_sd_alpha = 1.0; //standard deviation of the distribution of mutational effects for alpha
	mu_mean_beta = 0.0; //mean of the distribution of mutational effects for beta
	mu_sd_beta = 0.1; //standard deviation of the distribution of mutational effects for beta
	mu_mean_D0 = 0.0; //mean of the distribution of mutational effects for D0
	mu_sd_D0 = 0.1; //standard deviation of the distribution of mutational effects for D0
	//Parameters for density-dependent dispersal
	mu_mean_alphaD = 0.0; //mean of the distribution of mutational effects for alphaD
	mu_sd_alphaD = 1.0; //standard deviation of the distribution of mutational effects for alphaD
	mu_mean_betaD = 0.0; //mean of the distribution of mutational effects for betaD
	mu_sd_betaD = 0.1; //standard deviation of the distribution of mutational effects for betaD
	//Interaction between BC and density 
	mu_mean_gamma = 0.0; //mean of the distribution of mutational effects for gamma
	mu_sd_gamma = 0.1; //standard deviation of the distribution of mutational effects for gamma
	//Resource-dependent body condition
	alpha_res = 8.0; //slope of logistic function inflection point
	beta_res = 0.5; //inflection point
	//Reproduction
	max_fec = 8.0; //maximum fecundity for body-condition = 1
	BCfecundity = true; //BC-dependent fecundity
	//Dispersal
	disp_mort = 0.0; //constant dispersal mortality
	disp_function = 4;	// 0 = linear; 
						// 1 = body condition dep. (3-parameters logistic); 
						// 2 = body-condition and density dependent disp. (5 parameters)
						// 3 = density-dependent (3-parameters logistic); 
						// 4 = condition- (body or density) -independent dispersal
	disp_type = 0; // 0 = nearest_neighbour; 1 = global dispersal; 2 = kernel
	mean_dist = 2.0; //mean dispersal distance in the case of dispersal kernel
	//Survival dispersal cost
	disp_mort_type = 0; // 0 = constant (disp_mort); // 1 = BC-dependent (beta_surv, alpha_surv)
	beta_surv = 0.8; //inflection point BC-dependent dispersal survival
	alpha_surv = 0.0; //slope BC-dependent dispersal survival
	//Energetic dispersal cost
	disp_cost = 0.3; //dispersal cost paid as reduction in BC (hence fecundity for females)
	//Survival (after dispersal)
	survival = true; //yes/no
	//Environmental stochasticity
	ac = 0.2; //temporal autocorrelation (-9 = no env. stochasticity)
	std = 0.8; //standard deviation of env. fluctuations
	local_ext = 0.0; //local extinction probability
}


Parameters::~Parameters()
{
}
