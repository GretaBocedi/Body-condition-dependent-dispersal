#include "BCDispersal.h" 

//---------------------------------------------------------------------------
// Main function for running on the Linux Cluster
#if CLUSTER
int main(int argc, char* argv[])
{
	// Get the current directory.
	char* buffer = getcwd(NULL, 0);
	dir = buffer;
	free(buffer);
	dir = dir + "/"; //Current directory path
	dirOut = dir + "Outputs/"; //Outpus folder path

	RunModel();

	cout << "Simulation completed" << endl;

	return 0;
}
#else
int _tmain(int argc, _TCHAR* argv[]) // Main function for running on Windows
{
	// Get the current directory.
	char* buffer = _getcwd(NULL, 0);
	dir = buffer;
	free(buffer);
	dir = dir + "\\"; //Current directory path
	dirOut = dir + "Outputs\\"; //Outpus folder path

	RunModel();

	cout << "Simulation completed" << endl;

	return 0;
}
#endif
//---------------------------------------------------------------------------
// RunModel() handles looping through replicates, years and generations
void RunModel(void){
	cout << "Simulation nr. " << para.SimNr << endl;

	// create output files
	outPar();
	if (para.ind_int != 0) outInd_header();
	if (para.pop_int != 0) outPop_header();


	//Loop through REPLICATES --------------------------------------------
	for (r = 0; r < para.rep; r++) {
		cout << "rep = " << r << endl;

		std::normal_distribution<> normEnv(0.0, para.std);

		//Initialise --------------------------------------------------
		Initialisation(para.x_max, para.y_max);

		//Loop through GENERATIONS ------------------------------------
		for (g = 0; g < para.gen; g++) {

			if (para.ac > -9) env_stoch(normEnv); // environmental stochasticity 
			if (para.local_ext > 0.0) local_extinction(); //local extinction

			reproduction();
		
			if (totNoffs < 1) break;

			dispersal();

			if (para.survival) survival();

			//output populations
			if (g % para.pop_int == 0.0 || (g > para.expansion_start && g % 10 == 0)){
				for (int i = 0; i < para.x_max; i++) {
					for (int j = 0; j < para.y_max; j++)
						Land[i][j]->outPop();
				}
			}
			
			if (Nfemales < 1 || Nmales < 2) {
				cout << "Extinct" << endl;
				break;
			}
		} //END of the GENERATIONS loop

		//Delete landscape
		for (int i = 0; i < para.x_max; i++) {
			for (int j = 0; j < para.y_max; j++) {
				Land[i][j]->DeleteAdults();
				delete Land[i][j]; Land[i][j] = NULL;
			}
			delete[] Land[i]; Land[i] = NULL;
		}
		delete[] Land; Land = NULL;

	} //END of the REPLICATES loop

	// close output files
	if (para.ind_int != 0) inds.close();
	if (para.pop_int != 0) pops.close();
}
//---------------------------------------------------------------------------
// Converts integers into strings
const string Int2Str(const int x)
{
	ostringstream o;
	if (!(o << x)) return "ERROR";
	return o.str();
}
//---------------------------------------------------------------------------
// Converts doubles into strings
const string Float2Str(const double x)
{
	ostringstream o;
	if (!(o << x)) return "ERROR";
	return o.str();
}
//---------------------------------------------------------------------------
// Initialise the landscape, populations and individuals
void Initialisation(int xmax, int ymax){
	Individual *ind;
	vector<Individual>::iterator iter;

	Ntot = 0;
	Nfemales = 0;
	Nmales = 0;

	// uniform distribution to sample amount of resources in each cell
	std::uniform_real_distribution<> uniRes(para.minR, para.maxR);
	// normal distributions for allele initialisation
	std::normal_distribution<> distrAlpha(para.alpha_Imean / 2.0, para.alpha_Isd / sqrt(2.0));
	std::normal_distribution<> distrBeta(para.beta_Imean / 2.0, para.beta_Isd / sqrt(2.0));
	std::normal_distribution<> distrD0(para.D0_Imean / 2.0, para.D0_Isd / sqrt(2.0));
	std::normal_distribution<> distrAlphaD(para.alphaD_Imean / 2.0, para.alphaD_Isd / sqrt(2.0));
	std::normal_distribution<> distrBetaD(para.betaD_Imean / 2.0, para.betaD_Isd / sqrt(2.0));
	std::normal_distribution<> distrGamma(para.gamma_Imean / 2.0, para.gamma_Isd / sqrt(2.0));

	//Create landscape, assign resources and initialise individuals
	Land = new Cell **[xmax];
	for (int x = 0; x < xmax; x++) {
		Land[x] = new Cell *[ymax];
		for (int y = 0; y < ymax; y++) {
			Land[x][y] = new Cell(x, y);	
			Land[x][y]->Res = uniRes(rdgen); //assign a random amount of resources to the cell
			
			//initialise population within the cell
			if (y < para.max_seedY) {
				//Initialise a number of individuals (N) as if ideally everyone got 1 unit of resources
				Land[x][y]->N = (int)Land[x][y]->Res;

				//Initialise individuals
				for (int i = 0; i < Land[x][y]->N; i++) {
					ind = new Individual(sex(rdgen), x, y);

					ind->initialise(para.disp_function, distrAlpha, distrBeta, distrD0, distrAlphaD, distrBetaD, distrGamma);

					if (ind->sex) {
						Land[x][y]->females.push_back(*ind);
						Land[x][y]->Nfem++;
					}
					else {
						Land[x][y]->males.push_back(*ind);
						Land[x][y]->Nmale++;
					}
					delete ind;
				}

				if (Land[x][y]->N > 0 && Land[x][y]->Res > 0.0) {
					//Initial individuals acquire resources
					resource_competition(x, y, true);
				}

				//Update totals
				Ntot += Land[x][y]->N;
				Nfemales += Land[x][y]->Nfem;
				Nmales += Land[x][y]->Nmale;
			}
		}
	}

	maxY = para.max_seedY;
}
//---------------------------------------------------------------------------
// Applies environmental variation in the amount of resources in each cell of the landscape
void env_stoch(std::normal_distribution<> normEnv){
	for (int x = 0; x < para.x_max; x++){
		for (int y = 0; y < para.y_max; y++){
			Land[x][y]->eps = Land[x][y]->eps*para.ac + normEnv(rdgen)*sqrt(1.0 - para.ac*para.ac);
			Land[x][y]->Res = meanRes * (1.0 + Land[x][y]->eps);
			if (Land[x][y]->Res < para.minR) Land[x][y]->Res = para.minR;
			if (Land[x][y]->Res > para.maxR) Land[x][y]->Res = para.maxR;
		}
	}
}
//---------------------------------------------------------------------------
// For each cells in the landscape, determines whether the population goes extinct and, if yes, it deletes the individuals
void local_extinction(void){
	std::bernoulli_distribution ext(para.local_ext);

	for (int x = 0; x < para.x_max; x++) {
		for (int y = 0; y < maxY; y++) {
			if (Land[x][y]->N > 0) {
				if (ext(rdgen)) {
					Land[x][y]->DeleteAdults();
				}
			}
		}
	}
}
//---------------------------------------------------------------------------
// Individuals compete to acquire resources in the cell
// adult = true - competition among adults (initialisation)
// adult = fasle - competition among offspring 
void resource_competition(int x, int y, bool adult){
	int pars; // number of resource parcels
	double mean_n_par; //mean number of parcels per individual
	vector<Individual>::iterator iter;

	// Distribution of resource parcel size
	std::uniform_real_distribution<> uniParcel(para.minpar, para.maxpar);
	
	//number of expected parcels
	pars = (int)std::round(Land[x][y]->Res / ((para.minpar + para.maxpar) / 2.0));

	if (pars > 0.0){
		if (adult){
			mean_n_par = (double)pars / (double)Land[x][y]->N;

			//distribution to sample the nr. of parcels each individual acquires
			std::poisson_distribution<> parcels(mean_n_par);
			//loop through females
			for (iter = Land[x][y]->females.begin(); iter != Land[x][y]->females.end(); iter++){
				iter->n_par = parcels(rdgen);
				for (int i = 0; i < iter->n_par; i++) iter->resources += uniParcel(rdgen);
				iter->bodyCond(para.alpha_res, para.beta_res);
				if (iter->alive == false){
					Land[x][y]->N--;
					Land[x][y]->Nfem--;
				}
			}
			//loop through males
			for (iter = Land[x][y]->males.begin(); iter != Land[x][y]->males.end(); iter++){
				iter->n_par = parcels(rdgen);
				for (int i = 0; i < iter->n_par; i++) iter->resources += uniParcel(rdgen);
				iter->bodyCond(para.alpha_res, para.beta_res);
				if (iter->alive == false){
					Land[x][y]->N--;
					Land[x][y]->Nmale--;
				}
			}
		}
		else{ //offspring
			mean_n_par = (double)pars / (double)Land[x][y]->Noffs;

			//distribution to sample the nr. of parcels each individual acquires
			std::poisson_distribution<> parcels(mean_n_par);

			//loop through female offspring
			for (iter = Land[x][y]->Jfemales.begin(); iter != Land[x][y]->Jfemales.end(); iter++){
				iter->n_par = parcels(rdgen);
				for (int i = 0; i < iter->n_par; i++) iter->resources += uniParcel(rdgen);
				iter->bodyCond(para.alpha_res, para.beta_res);
				if (iter->alive) {
					Land[x][y]->survJfem.push_back(*iter);
					Land[x][y]->sum_BC += iter->BC;
				}
				else{
					iter->deleteInd();
					Land[x][y]->Noffs--;
					Land[x][y]->Foffs--;
				}

			}
			//loop through male offspring
			for (iter = Land[x][y]->Jmales.begin(); iter != Land[x][y]->Jmales.end(); iter++){
				iter->n_par = parcels(rdgen);
				for (int i = 0; i < iter->n_par; i++) iter->resources += uniParcel(rdgen);
				iter->bodyCond(para.alpha_res, para.beta_res);
				if (iter->alive) {
					Land[x][y]->survJmales.push_back(*iter);
					Land[x][y]->sum_BC += iter->BC;
				}
				else{
					iter->deleteInd();
					Land[x][y]->Noffs--;
					Land[x][y]->Moffs--;
				}
			}
		}
	}
	else{ //If there are offspring they all die
		if (adult == false) {
			//loop through female offspring
			for (iter = Land[x][y]->Jfemales.begin(); iter != Land[x][y]->Jfemales.end(); iter++) {
				iter->deleteInd();
				Land[x][y]->Noffs--;
				Land[x][y]->Foffs--;
			}
			//loop through male offspring
			for (iter = Land[x][y]->Jmales.begin(); iter != Land[x][y]->Jmales.end(); iter++) {
				iter->deleteInd();
				Land[x][y]->Noffs--;
				Land[x][y]->Moffs--;
			}
		}
	}
}
// ---------------------------------------------------------------------------
void reproduction(void){
	int dads;
	int dad; //index referring to positions in the males vector
	int offs;
	double fec;
	double sum; 
	double rdn;
	double *cum; //cumulative distribution of males probability of being chosen
	Individual *ind;
	vector<Individual>::iterator iter;

	std::uniform_real_distribution<> unif_real(0.0, 1.0);

	totNoffs = 0;
	totFoffs = 0;
	totMoffs = 0;

	//Loop through the landscape 
	for (int x = 0; x < para.x_max; x++){
		for (int y = 0; y < maxY; y++){

			//if there are individuals of both sexes in the cell proceed with reproduction
			if (Land[x][y]->Nfem > 0 && Land[x][y]->Nmale > 0){

				//Initialise cumulative distribution for sampling fathers (male probabliity of mating depends on their body-condition)
				dads = (int)Land[x][y]->males.size();
				cum = new double[dads];
				sum = 0.0;
				for (iter = Land[x][y]->males.begin(); iter != Land[x][y]->males.end(); iter++) sum += iter->BC;
				cum[0] = Land[x][y]->males[0].BC / sum;
				for (int i = 1; i < dads; i++) cum[i] = cum[i-1] + Land[x][y]->males[i].BC / sum;

				//Loop through females 
				for (iter = Land[x][y]->females.begin(); iter != Land[x][y]->females.end(); iter++){
					if (iter->alive){

						if (para.BCfecundity) fec = iter->BC * para.max_fec; //fecundity as a linear function of body-condition
						else fec = para.max_fec; //constant fecundity
				
						//sample the number of offspring (demographic stochasticity)
						std::poisson_distribution<> pois(fec);
						offs = pois(rdgen);

						//Initialise the offspring 
						for (int i = 0; i < offs; i++){
							//sample the father 
							//assume complete promiscuity, hence re-sample the father for each offspring
							dad = 0;
							do{
								rdn = unif_real(rdgen);
								for (int j = 0; j < dads; j++) {
									if (rdn <= cum[j]){
										dad = j;
										break;
									}
								}
							} while (Land[x][y]->males[dad].alive == false);

							//create new individual
							ind = new Individual(sex(rdgen), x, y);

							inheritance(ind, *iter, Land[x][y]->males[dad]); //inherit traits from the parents

							if (ind->sex){
								Land[x][y]->Jfemales.push_back(*ind);
								Land[x][y]->Foffs++;
							}
							else{
								Land[x][y]->Jmales.push_back(*ind);
								Land[x][y]->Moffs++;
							}
							delete ind;
							Land[x][y]->Noffs++;
						}
					}
				}

				delete[] cum; cum = NULL;

				if (Land[x][y]->Noffs > 0){
					//Offspring compete for resources
					resource_competition(x, y, false);
					Land[x][y]->Jfemales.clear();
					Land[x][y]->Jmales.clear();
				}
			}
			totNoffs += Land[x][y]->Noffs;
			totFoffs += Land[x][y]->Foffs;
			totMoffs += Land[x][y]->Moffs;

			//Clear adults
			Land[x][y]->DeleteAdults();
		}
	}

	Ntot = 0;
	Nfemales = 0;
	Nmales = 0;

	// Mutations
	if (para.mu > 0.0 && totNoffs > 0) mutations();
}
//---------------------------------------------------------------------------
void inheritance(Individual *pup, Individual mom, Individual dad){
	int rdn, rdn2;

	switch (para.disp_function){
	case 0: //linear function
		//alpha
		pup->traits.alpha = new double[2];
		pup->traits.g_alpha = new double(0.0);
		pup->traits.p_alpha = new double(0.0);
		//beta
		pup->traits.beta = new double[2];
		pup->traits.g_beta = new double(0.0);
		pup->traits.p_beta = new double(0.0);
		//Inherit alleles from mather and father
		rdn = bern(rdgen);
		rdn2 = bern(rdgen);
		pup->traits.alpha[0] = mom.traits.alpha[rdn];
		pup->traits.alpha[1] = dad.traits.alpha[rdn2];

		rdn = bern(rdgen);
		rdn2 = bern(rdgen);
		pup->traits.beta[0] = mom.traits.beta[rdn];
		pup->traits.beta[1] = dad.traits.beta[rdn2];
		//Genotypes (sum of alleles)
		*pup->traits.g_alpha = pup->traits.alpha[0] + pup->traits.alpha[1];
		*pup->traits.g_beta = pup->traits.beta[0] + pup->traits.beta[1];
		//Phenotypes
		*pup->traits.p_alpha = *pup->traits.g_alpha;
		*pup->traits.p_beta = *pup->traits.g_beta;
		break;
	case 1: //body condition dep. (3-parameters logistic); 
		//alpha
		pup->traits.alpha = new double[2];
		pup->traits.g_alpha = new double(0.0);
		pup->traits.p_alpha = new double(0.0);
		//beta
		pup->traits.beta = new double[2];
		pup->traits.g_beta = new double(0.0);
		pup->traits.p_beta = new double(0.0);
		//D0
		pup->traits.D0 = new double[2];
		pup->traits.g_D0 = new double(0.0);
		pup->traits.p_D0 = new double(0.0);
		//Inherit alleles from mather and father
		rdn = bern(rdgen);
		rdn2 = bern(rdgen);
		pup->traits.alpha[0] = mom.traits.alpha[rdn];
		pup->traits.alpha[1] = dad.traits.alpha[rdn2];

		rdn = bern(rdgen);
		rdn2 = bern(rdgen);
		pup->traits.beta[0] = mom.traits.beta[rdn];
		pup->traits.beta[1] = dad.traits.beta[rdn2];

		rdn = bern(rdgen);
		rdn2 = bern(rdgen);
		pup->traits.D0[0] = mom.traits.D0[rdn];
		pup->traits.D0[1] = dad.traits.D0[rdn2];
		//Genotypes (sum of alleles)
		*pup->traits.g_alpha = pup->traits.alpha[0] + pup->traits.alpha[1];
		*pup->traits.g_beta = pup->traits.beta[0] + pup->traits.beta[1];
		*pup->traits.g_D0 = pup->traits.D0[0] + pup->traits.D0[1];
		//Phenotypes
		*pup->traits.p_alpha = *pup->traits.g_alpha;
		*pup->traits.p_beta = *pup->traits.g_beta;
		*pup->traits.p_D0 = *pup->traits.g_D0;
		if (*pup->traits.p_D0 < 0.0) *pup->traits.p_D0 = 0.0;// Constrain D0 between 0 and 1
		if (*pup->traits.p_D0 > 1.0) *pup->traits.p_D0 = 1.0;
		break;
	case 2: //body condition and density dependent disp. (5 parameters)
		//alpha
		pup->traits.alpha = new double[2];
		pup->traits.g_alpha = new double(0.0);
		pup->traits.p_alpha = new double(0.0);
		//beta
		pup->traits.beta = new double[2];
		pup->traits.g_beta = new double(0.0);
		pup->traits.p_beta = new double(0.0);
		//alphaD
		pup->traits.alphaD = new double[2];
		pup->traits.g_alphaD = new double(0.0);
		pup->traits.p_alphaD = new double(0.0);
		//betaD
		pup->traits.betaD = new double[2];
		pup->traits.g_betaD = new double(0.0);
		pup->traits.p_betaD = new double(0.0);
		//gamma
		pup->traits.gamma = new double[2];
		pup->traits.g_gamma = new double(0.0);
		pup->traits.p_gamma = new double(0.0);
		//Inherit alleles from mather and father
		rdn = bern(rdgen);
		rdn2 = bern(rdgen);
		pup->traits.alpha[0] = mom.traits.alpha[rdn];
		pup->traits.alpha[1] = dad.traits.alpha[rdn2];

		rdn = bern(rdgen);
		rdn2 = bern(rdgen);
		pup->traits.beta[0] = mom.traits.beta[rdn];
		pup->traits.beta[1] = dad.traits.beta[rdn2];

		rdn = bern(rdgen);
		rdn2 = bern(rdgen);
		pup->traits.alphaD[0] = mom.traits.alphaD[rdn];
		pup->traits.alphaD[1] = dad.traits.alphaD[rdn2];

		rdn = bern(rdgen);
		rdn2 = bern(rdgen);
		pup->traits.betaD[0] = mom.traits.betaD[rdn];
		pup->traits.betaD[1] = dad.traits.betaD[rdn2];

		rdn = bern(rdgen);
		rdn2 = bern(rdgen);
		pup->traits.gamma[0] = mom.traits.gamma[rdn];
		pup->traits.gamma[1] = dad.traits.gamma[rdn2];
		//Genotypes (sum of alleles)
		*pup->traits.g_alpha = pup->traits.alpha[0] + pup->traits.alpha[1];
		*pup->traits.g_beta = pup->traits.beta[0] + pup->traits.beta[1];
		*pup->traits.g_alphaD = pup->traits.alphaD[0] + pup->traits.alphaD[1];
		*pup->traits.g_betaD = pup->traits.betaD[0] + pup->traits.betaD[1];
		*pup->traits.g_gamma = pup->traits.gamma[0] + pup->traits.gamma[1];
		//Phenotypes
		*pup->traits.p_alpha = *pup->traits.g_alpha;
		*pup->traits.p_beta = *pup->traits.g_beta;
		*pup->traits.p_alphaD = *pup->traits.g_alphaD;
		*pup->traits.p_betaD = *pup->traits.g_betaD;
		*pup->traits.p_gamma = *pup->traits.g_gamma;
		break;
	case 3: //density-dependent disp. (3-parameters logistic); 
		//alphaD
		pup->traits.alphaD = new double[2];
		pup->traits.g_alphaD = new double(0.0);
		pup->traits.p_alphaD = new double(0.0);
		//betaD
		pup->traits.betaD = new double[2];
		pup->traits.g_betaD = new double(0.0);
		pup->traits.p_betaD = new double(0.0);
		//D0
		pup->traits.D0 = new double[2];
		pup->traits.g_D0 = new double(0.0);
		pup->traits.p_D0 = new double(0.0);
		//Inherit alleles from mather and father
		rdn = bern(rdgen);
		rdn2 = bern(rdgen);
		pup->traits.alphaD[0] = mom.traits.alphaD[rdn];
		pup->traits.alphaD[1] = dad.traits.alphaD[rdn2];

		rdn = bern(rdgen);
		rdn2 = bern(rdgen);
		pup->traits.betaD[0] = mom.traits.betaD[rdn];
		pup->traits.betaD[1] = dad.traits.betaD[rdn2];

		rdn = bern(rdgen);
		rdn2 = bern(rdgen);
		pup->traits.D0[0] = mom.traits.D0[rdn];
		pup->traits.D0[1] = dad.traits.D0[rdn2];
		//Genotypes (sum of alleles)
		*pup->traits.g_alphaD = pup->traits.alphaD[0] + pup->traits.alphaD[1];
		*pup->traits.g_betaD = pup->traits.betaD[0] + pup->traits.betaD[1];
		*pup->traits.g_D0 = pup->traits.D0[0] + pup->traits.D0[1];
		//Phenotypes
		*pup->traits.p_alphaD = *pup->traits.g_alphaD;
		*pup->traits.p_betaD = *pup->traits.g_betaD;
		*pup->traits.p_D0 = *pup->traits.g_D0;
		if (*pup->traits.p_D0 < 0.0) *pup->traits.p_D0 = 0.0;// Constrain D0 between 0 and 1
		if (*pup->traits.p_D0 > 1.0) *pup->traits.p_D0 = 1.0;
		break;
	case 4: //condition independent dispersal
		//D0
		pup->traits.D0 = new double[2];
		pup->traits.g_D0 = new double(0.0);
		pup->traits.p_D0 = new double(0.0);
		//Inherit alleles from mather and father
		rdn = bern(rdgen);
		rdn2 = bern(rdgen);
		pup->traits.D0[0] = mom.traits.D0[rdn];
		pup->traits.D0[1] = dad.traits.D0[rdn2];
		//Genotypes (sum of alleles)
		*pup->traits.g_D0 = pup->traits.D0[0] + pup->traits.D0[1];
		//Phenotypes
		*pup->traits.p_D0 = *pup->traits.g_D0;
		if (*pup->traits.p_D0 < 0.0) *pup->traits.p_D0 = 0.0;// Constrain D0 between 0 and 1
		if (*pup->traits.p_D0 > 1.0) *pup->traits.p_D0 = 1.0;
		break;
	}

}
//---------------------------------------------------------------------------
void mutations(void){
	bool vivo;
	int Nmu, ind, allele, x, y;

	std::poisson_distribution<> pois(2 * para.mu * (double)totNoffs); //distribution to sample the total nr. of mutations in the population

	if (para.disp_function < 3){
		//Mutations in Alpha
		Nmu = pois(rdgen); //total nr. of mutations
		for (int i = 0; i < Nmu; i++) {
			do{
				x = sample_x(rdgen); //sample cell's coordinates
				y = sample_y(rdgen);
			} while (Land[x][y]->Noffs < 1);
			std::uniform_int_distribution<> unif(0, Land[x][y]->Noffs - 1);
			do{
				ind = unif(rdgen); //sample individual index (survJfem + survJmales)
				if (ind < Land[x][y]->Foffs) vivo = Land[x][y]->survJfem[ind].alive;
				else vivo = Land[x][y]->survJmales[ind - Land[x][y]->Foffs].alive;
			} while (vivo == false);
			allele = bern(rdgen); //sample allele
			if (ind < Land[x][y]->Foffs) Land[x][y]->survJfem[ind].mutation(0, allele);
			else Land[x][y]->survJmales[ind - Land[x][y]->Foffs].mutation(0, allele);
		}
		//Mutations in Beta
		Nmu = pois(rdgen); //total nr. of mutations
		for (int i = 0; i < Nmu; i++) {
			do{
				x = sample_x(rdgen); //sample cell's coordinates
				y = sample_y(rdgen);
			} while (Land[x][y]->Noffs < 1);
			std::uniform_int_distribution<> unif(0, Land[x][y]->Noffs - 1);
			do{
				ind = unif(rdgen); //sample individual index (survJfem + survJmales)
				if (ind < Land[x][y]->Foffs) vivo = Land[x][y]->survJfem[ind].alive;
				else vivo = Land[x][y]->survJmales[ind - Land[x][y]->Foffs].alive;
			} while (vivo == false);
			allele = bern(rdgen); //sample allele
			if (ind < Land[x][y]->Foffs) Land[x][y]->survJfem[ind].mutation(1, allele);
			else Land[x][y]->survJmales[ind - Land[x][y]->Foffs].mutation(1, allele);
		}
	}

	if (para.disp_function == 1 || para.disp_function > 2){
		//Mutations in D0
			Nmu = pois(rdgen); //total nr. of mutations
			for (int i = 0; i < Nmu; i++) {
				do{
					x = sample_x(rdgen); //sample cell's coordinates
					y = sample_y(rdgen);
				} while (Land[x][y]->Noffs < 1);
				std::uniform_int_distribution<> unif(0, Land[x][y]->Noffs - 1);
				do{
					ind = unif(rdgen); //sample individual index (survJfem + survJmales)
					if (ind < Land[x][y]->Foffs) vivo = Land[x][y]->survJfem[ind].alive;
					else vivo = Land[x][y]->survJmales[ind - Land[x][y]->Foffs].alive;
				} while (vivo == false);
				allele = bern(rdgen); //sample allele
				if (ind < Land[x][y]->Foffs) Land[x][y]->survJfem[ind].mutation(2, allele);
				else Land[x][y]->survJmales[ind - Land[x][y]->Foffs].mutation(2, allele);
			}
	}

	if (para.disp_function == 2 || para.disp_function == 3){
		//Mutations in AlphaD
		Nmu = pois(rdgen); //total nr. of mutations
		for (int i = 0; i < Nmu; i++) {
			do{
				x = sample_x(rdgen); //sample cell's coordinates
				y = sample_y(rdgen);
			} while (Land[x][y]->Noffs < 1);
			std::uniform_int_distribution<> unif(0, Land[x][y]->Noffs - 1);
			do{
				ind = unif(rdgen); //sample individual index (survJfem + survJmales)
				if (ind < Land[x][y]->Foffs) vivo = Land[x][y]->survJfem[ind].alive;
				else vivo = Land[x][y]->survJmales[ind - Land[x][y]->Foffs].alive;
			} while (vivo == false);
			allele = bern(rdgen); //sample allele
			if (ind < Land[x][y]->Foffs) Land[x][y]->survJfem[ind].mutation(3, allele);
			else Land[x][y]->survJmales[ind - Land[x][y]->Foffs].mutation(3, allele);
		}

		//Mutations in BetaD
		Nmu = pois(rdgen); //total nr. of mutations
		for (int i = 0; i < Nmu; i++) {
			do{
				x = sample_x(rdgen); //sample cell's coordinates
				y = sample_y(rdgen);
			} while (Land[x][y]->Noffs < 1);
			std::uniform_int_distribution<> unif(0, Land[x][y]->Noffs - 1);
			do{
				ind = unif(rdgen); //sample individual index (survJfem + survJmales)
				if (ind < Land[x][y]->Foffs) vivo = Land[x][y]->survJfem[ind].alive;
				else vivo = Land[x][y]->survJmales[ind - Land[x][y]->Foffs].alive;
			} while (vivo == false);
			allele = bern(rdgen); //sample allele
			if (ind < Land[x][y]->Foffs) Land[x][y]->survJfem[ind].mutation(4, allele);
			else Land[x][y]->survJmales[ind - Land[x][y]->Foffs].mutation(4, allele);
		}
	}
	if (para.disp_function == 2){
		//Mutations in Gamma
		Nmu = pois(rdgen); //total nr. of mutations
		for (int i = 0; i < Nmu; i++) {
			do{
				x = sample_x(rdgen); //sample cell's coordinates
				y = sample_y(rdgen);
			} while (Land[x][y]->Noffs < 1);
			std::uniform_int_distribution<> unif(0, Land[x][y]->Noffs - 1);
			do{
				ind = unif(rdgen); //sample individual index (survJfem + survJmales)
				if (ind < Land[x][y]->Foffs) vivo = Land[x][y]->survJfem[ind].alive;
				else vivo = Land[x][y]->survJmales[ind - Land[x][y]->Foffs].alive;
			} while (vivo == false);
			allele = bern(rdgen); //sample allele
			if (ind < Land[x][y]->Foffs) Land[x][y]->survJfem[ind].mutation(5, allele);
			else Land[x][y]->survJmales[ind - Land[x][y]->Foffs].mutation(5, allele);
		}
	}
}
//---------------------------------------------------------------------------
// Function declared in the class Individual. Applies a specific mutation to the individual
void Individual::mutation(int trait, int allele){

	switch (trait){
	case 0: //alpha
		*traits.g_alpha -= traits.alpha[allele];
		traits.alpha[allele] += normAlpha(rdgen);
		*traits.g_alpha += traits.alpha[allele];
		*traits.p_alpha = *traits.g_alpha;
		break;
	case 1: //Beta
		*traits.g_beta -= traits.beta[allele];
		traits.beta[allele] += normBeta(rdgen);
		*traits.g_beta += traits.beta[allele];
		*traits.p_beta = *traits.g_beta;
		break;
	case 2: //D0
		*traits.g_D0 -= traits.D0[allele];
		traits.D0[allele] += normD0(rdgen);
		*traits.g_D0 += traits.D0[allele];
		*traits.p_D0 = *traits.g_D0;
		if (*traits.p_D0 < 0.0) *traits.p_D0 = 0.0;
		if (*traits.p_D0 > 1.0) *traits.p_D0 = 1.0;
		break;
	case 3: //alphaD
		*traits.g_alphaD -= traits.alphaD[allele];
		traits.alphaD[allele] += normAlphaD(rdgen);
		*traits.g_alphaD += traits.alphaD[allele];
		*traits.p_alphaD = *traits.g_alphaD;
		break;
	case 4: //BetaD
		*traits.g_betaD -= traits.betaD[allele];
		traits.betaD[allele] += normBetaD(rdgen);
		*traits.g_betaD += traits.betaD[allele];
		*traits.p_betaD = *traits.g_betaD;
		break;
	case 5: //gamma
		*traits.g_gamma -= traits.gamma[allele];
		traits.gamma[allele] += normGamma(rdgen);
		*traits.g_gamma += traits.gamma[allele];
		*traits.p_gamma = *traits.g_gamma;
		break;
	}
}
//---------------------------------------------------------------------------
void dispersal(void){
	vector<Individual>::iterator iter;
	double surv;

	int prev_maxY = maxY; 

	//Loop through the landscape 
	for (int x = 0; x < para.x_max; x++){
		for (int y = 0; y < prev_maxY; y++){

			Land[x][y]->emigrants = 0;

			if (Land[x][y]->Noffs > 0){

				//Loop through female offspring
				for (iter = Land[x][y]->survJfem.begin(); iter != Land[x][y]->survJfem.end(); iter++){
					if (iter->alive){
						//determine the individual dispersal probability
						iter->dispersalP(para.disp_function, (double)Land[x][y]->Noffs / Land[x][y]->Res);
						if (iter->d < 0.000000000001) iter->d = 0.000000000001;

						//determine whether it disperses
						std::bernoulli_distribution disp(iter->d);
						if (disp(rdgen)) {
							iter->disperse(para.x_max, para.y_max, para.disp_type);
														
							Land[x][y]->emigrants++;
							Land[x][y]->sum_BC -= iter->BC;
							
							if (para.disp_mort_type == 0){ //Constant dispersal mortality 
								surv = 1.0 - para.disp_mort;
							}
							else{ // Body condition-dependent dispersal survival probability
								surv = 1.0 / (1.0 + exp(-(iter->BC - para.beta_surv)*para.alpha_surv));								
							}

							std::bernoulli_distribution dispSurv(surv);
							if (dispSurv(rdgen)) { //dispersal mortality
								//energetic cost of dispersal
								if (para.disp_mort_type == 0 && para.disp_cost > 0.0){
									iter->BC -= para.disp_cost; //dispersal cost to BC
									//if BC <= 0 the individual dies
									if (iter->BC <= 0.0) iter->alive = false;
								}
							}
							else iter->alive = false;

							if (iter->alive == true){
								//add the individual to the adult vector in its current cell							
								Land[iter->new_x][iter->new_y]->females.push_back(*iter);
								Land[iter->new_x][iter->new_y]->Nfem++;
								Land[iter->new_x][iter->new_y]->N++;
								Land[iter->new_x][iter->new_y]->sum_BC += iter->BC;
								Ntot++;
								Nfemales++;
								if (g % para.ind_int == 0.0) iter->outInd();

								if (iter->new_y + 1 > maxY) maxY = iter->new_y + 1;
							}
							else {
								if (g % para.ind_int == 0.0) iter->outInd();
								iter->deleteInd();
							}
						}
						else {
							Land[iter->x][iter->y]->females.push_back(*iter);
							Land[iter->x][iter->y]->Nfem++;
							Land[iter->x][iter->y]->N++;
							Ntot++;
							Nfemales++;
							if (g % para.ind_int == 0.0) iter->outInd();
						}
						
					}
					else iter->deleteInd();
				}
				//Loop through male offspring
				for (iter = Land[x][y]->survJmales.begin(); iter != Land[x][y]->survJmales.end(); iter++){
					if (iter->alive){
						//determine the individual dispersal probability
						iter->dispersalP(para.disp_function, (double)Land[x][y]->Noffs / Land[x][y]->Res);
						if (iter->d < 0.000000000001) iter->d = 0.000000000001;

						//determine whether it disperses
						std::bernoulli_distribution disp(iter->d);
						if (disp(rdgen)) {
							iter->disperse(para.x_max, para.y_max, para.disp_type);

							Land[x][y]->emigrants++;
							Land[x][y]->sum_BC -= iter->BC;

							if (para.disp_mort_type == 0){ //Constant dispersal mortality 
								surv = 1.0 - para.disp_mort;
							}
							else{ // Body condition-dependent dispersal survival probability
								surv = 1.0 / (1.0 + exp(-(iter->BC - para.beta_surv)*para.alpha_surv));
							}

							std::bernoulli_distribution dispSurv(surv);
							if (dispSurv(rdgen)) { //dispersal mortality
								//energetic cost of dispersal
								if (para.disp_mort_type == 0 && para.disp_cost > 0.0){
									iter->BC -= para.disp_cost; //dispersal cost to BC
									//if BC <= 0 the individual dies
									if (iter->BC <= 0.0) iter->alive = false;
								}
							}
							else iter->alive = false;

							if (iter->alive == true){
								//add the individual to the adult vector in its current cell							
								Land[iter->new_x][iter->new_y]->males.push_back(*iter);
								Land[iter->new_x][iter->new_y]->Nmale++;
								Land[iter->new_x][iter->new_y]->N++;
								Land[iter->new_x][iter->new_y]->sum_BC += iter->BC;
								Ntot++;
								Nmales++;
								if (g % para.ind_int == 0.0) iter->outInd();
								if (iter->new_y + 1 > maxY) maxY = iter->new_y + 1;
							}
							else {
								if (g % para.ind_int == 0.0) iter->outInd();
								iter->deleteInd();
							}
						}
						else {
							Land[iter->x][iter->y]->males.push_back(*iter);
							Land[iter->x][iter->y]->Nmale++;
							Land[iter->x][iter->y]->N++;
							Ntot++;
							Nmales++;
							if (g % para.ind_int == 0.0) iter->outInd();
						}
					}
					else iter->deleteInd();
				}
				Land[x][y]->survJfem.clear();
				Land[x][y]->survJmales.clear();

				Land[x][y]->prop_emigrants = (double)Land[x][y]->emigrants / (double)Land[x][y]->Noffs;

				Land[x][y]->Noffs = 0;
				Land[x][y]->Foffs = 0;
				Land[x][y]->Moffs = 0;
			}
		}
	}
	totNoffs = 0;
	totMoffs = 0;
	totFoffs = 0;
}
//---------------------------------------------------------------------------
// Function declared in the class Individual. Determines the new location of a dispersing individual 
void Individual::disperse(int max_x, int max_y, int type) {
	double x_rand, y_rand;
	double R1, dist, rndAngle;

	if (type == 0) {//nearest neighbour dispersal
		std::uniform_int_distribution<> sample_x(x - 1, x + 1);
		std::uniform_int_distribution<> sample_y(y - 1, y + 1);
		do {
			new_x = sample_x(rdgen);
			new_y = sample_y(rdgen);
		} while ((new_x == x && new_y == y) || new_x < 0 || new_x >(max_x - 1) || new_y < 0 || new_y >(max_y - 1));
	}
	else {
		if (type == 1) {//global dispersal
			std::uniform_int_distribution<> sample_x2(0, max_x - 1);
			std::uniform_int_distribution<> sample_y2(0, max_y - 1);
			do {
				new_x = sample_x2(rdgen);
				new_y = sample_y2(rdgen);
			} while ((new_x == x && new_y == y) || new_x < 0 || new_x >(max_x - 1) || new_y < 0 || new_y >(max_y - 1));
		}
		else {//negative-exponential kernel
			std::uniform_real_distribution<> unireal_disp(0.0, 0.999);
			std::uniform_real_distribution<> unireal_dispB(0.0000001, 1.0);
			//sample new location
			x_rand = unireal_disp(rdgen);
			y_rand = unireal_disp(rdgen);
			do {
				do {
					R1 = unireal_dispB(rdgen);
					dist = (-1.0 * para.mean_dist) * std::log(R1);
					rndAngle = unireal(rdgen) * 2.0 * PI;
					new_x = (int)(dist * cos(rndAngle) + x_rand + x);
					new_y = (int)(dist * sin(rndAngle) + y_rand + y);
				} while (new_x == x && new_y == y);
			} while (new_x < 0 || new_x >(max_x - 1) || new_y < 0 || new_y >(max_y - 1));
		}
	}

	dispersed = true;
}
//---------------------------------------------------------------------------
//BC-dependent survival (post dispersal)
void survival(void) {
	double surv;
	vector<Individual>::iterator iter; 

	//Loop through the landscape 
	for (int x = 0; x < para.x_max; x++) {
		for (int y = 0; y < maxY; y++) {

			if (Land[x][y]->N > 0) {

				//Loop through females
				for (iter = Land[x][y]->females.begin(); iter != Land[x][y]->females.end(); iter++) {
					if (iter->alive) {
						
						//determine survival probability based on density and condition
						surv = std::fmin(Land[x][y]->Res*(iter->BC / Land[x][y]->sum_BC),1.0);
						//determine whether it survives
						std::bernoulli_distribution survive(surv);
						if (survive(rdgen) == false) {
							iter->alive = false;
							Land[x][y]->Nfem--;
							Land[x][y]->N--;
						}
					}
				}
				//Loop through males
				for (iter = Land[x][y]->males.begin(); iter != Land[x][y]->males.end(); iter++) {
					if (iter->alive) {

						//determine survival probability based on density and condition
						surv = std::fmin(Land[x][y]->Res*(iter->BC / Land[x][y]->sum_BC), 1.0);
						//determine whether it survives
						std::bernoulli_distribution survive(surv);
						if (survive(rdgen) == false) {
							iter->alive = false;
							Land[x][y]->Nmale--;
							Land[x][y]->N--;
						}
					}
				}
			}
			Land[x][y]->sum_BC = 0.0;
		}
	}
}

// Output parameters --------------------------------------------------------
void outPar(void){
	string name;

	name = dirOut + "Sim" + Int2Str(para.SimNr) + "_Parameters.txt";
	par.open(name.c_str());
	par << "rep\t" << para.rep << endl;
	par << "gen\t" << para.gen << endl;
	par << "x_max\t" << para.x_max << endl << "y_max\t" << para.y_max << endl;
	par << "x_seed_max\t" << para.max_seedX << endl;
	par << "y_seed_max\t" << para.max_seedY << endl;
	par << "minResources\t" << para.minR<< endl;
	par << "maxResources\t" << para.maxR << endl;
	par << "minResParcelSize\t" << para.minpar << endl;
	par << "maxResParcelSize\t" << para.maxpar << endl;
	par << "max_fecundity\t" << para.max_fec<< endl;
	par << "BC-dependent fecundity\t" << para.BCfecundity << endl;
	par << "alpha_Imean\t" << para.alpha_Imean << endl;
	par << "alpha_Isd\t" << para.alpha_Isd << endl;
	par << "beta_Imean\t" << para.beta_Imean << endl;
	par << "beta_Isd\t" << para.beta_Isd << endl;
	par << "D0_Imean\t" << para.D0_Imean << endl;
	par << "D0_Isd\t" << para.D0_Isd << endl;
	par << "alphaD_Imean\t" << para.alphaD_Imean << endl;
	par << "alphaD_Isd\t" << para.alphaD_Isd << endl;
	par << "betaD_Imean\t" << para.betaD_Imean << endl;
	par << "betaD_Isd\t" << para.betaD_Isd << endl;
	par << "gamma_Imean\t" << para.gamma_Imean << endl;
	par << "gamma_Isd\t" << para.gamma_Isd << endl;
	par << "mutation probability\t" << para.mu << endl;
	par << "mu_mean_alpha\t" << para.mu_mean_alpha << endl;
	par << "mu_sd_alpha\t" << para.mu_sd_alpha << endl;
	par << "mu_mean_beta\t" << para.mu_mean_beta << endl;
	par << "mu_sd_beta\t" << para.mu_sd_beta << endl;
	par << "mu_mean_D0\t" << para.mu_mean_D0 << endl;
	par << "mu_sd_D0\t" << para.mu_sd_D0 << endl;
	par << "mu_mean_alphaD\t" << para.mu_mean_alpha << endl;
	par << "mu_sd_alphaD\t" << para.mu_sd_alpha << endl;
	par << "mu_mean_betaD\t" << para.mu_mean_beta << endl;
	par << "mu_sd_betaD\t" << para.mu_sd_beta << endl;
	par << "mu_mean_gamma\t" << para.mu_mean_D0 << endl;
	par << "mu_sd_gamma\t" << para.mu_sd_D0 << endl;
	par << "alpha_res\t" << para.alpha_res << endl;
	par << "beta_res\t" << para.beta_res << endl;
	par << "Dispersal function\t" << para.disp_function << endl;
	par << "Dispersal mortality type\t" << para.disp_mort_type << endl;
	if (para.disp_mort_type == 0) {
		par << "Dispersal mortality\t" << para.disp_mort << endl;
		par << "Dispersal energetic cost\t" << para.disp_cost << endl;
	}
	else par << "beta_surv\t" << para.beta_surv << endl << "alpha_surv\t" << para.alpha_surv << endl;
	par << "Transfer\t" << para.disp_type << endl;
	if (para.disp_type == 2) par << "mean dispersal distance\t" << para.mean_dist << endl;
	par << "Mortality after dispersal\t" << para.survival << endl;
	par << "Temporal autocorrelation\t" << para.ac << endl;
	par << "Temporal standard deviation\t" << para.std << endl;
	par << "Local extinction probability\t" << para.local_ext << endl;

	par.close();
}
//------------------------------------------------------------------------------
// Creates the Individuals output file and write the column headers
void outInd_header(void){
	string name = dirOut + "Sim" + Int2Str(para.SimNr) + "_Inds.txt";
	inds.open(name.c_str());

	//inds << "rep\tgen\tx\ty\tsex\tres\tn_par\tBC";
	//if (para.disp_function < 3)	inds << "\tp_beta\tp_alpha";
	//if (para.disp_function == 1 || para.disp_function > 2) inds << "\tp_D0";
	//if (para.disp_function == 2) inds << "\tp_betaD\tp_alphaD\tgamma";
	//if (para.disp_function == 3) inds << "\tp_betaD\tp_alphaD";
	//inds << "\td\tdispersed\talive\tnew_x\tnew_y\tdensNatal" << endl;

	inds << "rep\tgen\tx\ty\tBC";
	if (para.disp_function < 3)	inds << "\tp_beta\tp_alpha";
	if (para.disp_function == 1 || para.disp_function > 2) inds << "\tp_D0";
	if (para.disp_function == 2) inds << "\tp_betaD\tp_alphaD\tgamma";
	if (para.disp_function == 3) inds << "\tp_betaD\tp_alphaD";
	inds << "\td\tdispersed\talive\tdensNatal" << endl;
}
//------------------------------------------------------------------------------
// Function declared in the class Individual. Writes the outputs for each individual
void Individual::outInd(void){

	/*inds << r << "\t" << g << "\t" << x << "\t" << y << "\t" << sex << "\t" << resources << "\t" << n_par << "\t" << BC << "\t";
	if (para.disp_function < 3) inds << *traits.p_beta << "\t" << *traits.p_alpha << "\t";
	if(para.disp_function == 1 || para.disp_function > 2) inds << *traits.p_D0 << "\t";
	if (para.disp_function == 2) inds << *traits.p_betaD << "\t" << *traits.p_alphaD << "\t" << *traits.p_gamma << "\t";
	if (para.disp_function == 3) inds << *traits.p_betaD << "\t" << *traits.p_alphaD << "\t";
	inds << d << "\t" << dispersed << "\t" << alive << "\t" << new_x << "\t" << new_y << "\t" << (double)Land[x][y]->Noffs / Land[x][y]->Res << endl;*/

	inds << r << "\t" << g;
	if (dispersed) inds << "\t" << new_x << "\t" << new_y;
	else inds << "\t" << x << "\t" << y;
	inds << "\t" << BC << "\t";
	if (para.disp_function < 3) inds << *traits.p_beta << "\t" << *traits.p_alpha << "\t";
	if (para.disp_function == 1 || para.disp_function > 2) inds << *traits.p_D0 << "\t";
	if (para.disp_function == 2) inds << *traits.p_betaD << "\t" << *traits.p_alphaD << "\t" << *traits.p_gamma << "\t";
	if (para.disp_function == 3) inds << *traits.p_betaD << "\t" << *traits.p_alphaD << "\t";
	inds << d << "\t" << dispersed << "\t" << alive << "\t" << (double)Land[x][y]->Noffs / Land[x][y]->Res << endl;

}
//------------------------------------------------------------------------------
// Creates the Populations output file and write the column headers
void outPop_header(void){
	string name = dirOut + "Sim" + Int2Str(para.SimNr) + "_Pops.txt";
	pops.open(name.c_str());

	pops << "rep\tgen\tx\ty\tRes\tN\temigrants\tprop_emig" << endl;
}
//------------------------------------------------------------------------------
// Function declared in the class Cell. Writes the outputs for each population
void Cell::outPop(void){

	pops << r << "\t" << g << "\t" << x << "\t" << y << "\t" << Res << "\t" << N << "\t" << emigrants << "\t" << prop_emigrants << endl;
}
