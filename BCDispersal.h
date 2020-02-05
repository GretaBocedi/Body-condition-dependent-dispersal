/*------------------------------------------------------------------------------
Body condition-dependent dispersal: Model

Implements the main functions of the model.
Generates random number distributions.

THE JOINT EVOLUTION OF DENSITY- AND BODY CONDITION-DEPENDENT DISPERSAL
Baines, C.B., Travis, J.M.J, McCauley, S.J. & Bocedi, G. - submitted

Author: Greta Bocedi

Last updated: 29/01/2020 by G. Bocedi
--------------------------------------------------------------------------------*/

#pragma once


#define CLUSTER 0 //set to 1 to run on a Linux Cluster

#include <stdio.h>
#include <stdlib.h>
#if CLUSTER 
#include <unistd.h>
#else
#include <tchar.h> 
#include <direct.h>
#include <io.h>
#endif
#include <vector>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <math.h>
#include <numeric>
#include <time.h>
#include <random>
#include <iterator>

#include "Parameters.h"
#include "Individual.h"
#include "Cell.h"

using namespace std;


const double PI = 3.141592654;

//Global Objects and variables
Parameters para;
Cell ***Land;

string dir, dirOut;
std::ofstream par, inds, pops;
int maxY; //to track the maximum occupied y coordinate at each generation
int r, g; // counters for replicates and generations
int Ntot, Nfemales, Nmales, totNoffs, totFoffs, totMoffs; //totals

double meanRes = (para.minR + para.maxR) / 2.0; //mean resources

// Random numbers generators
// seed random number generator
std::random_device rd;
std::mt19937 rdgen(rd());

std::bernoulli_distribution sex(0.5);
std::bernoulli_distribution bern(0.5);
std::uniform_int_distribution<> sample_x(0, para.x_max - 1);
std::uniform_int_distribution<> sample_y(0, para.y_max - 1);
std::uniform_real_distribution<> unireal(0.0, 1.0);
std::bernoulli_distribution Mutate(para.mu); // Mutation probability
// normal distributions for mutational effects
std::normal_distribution<> normAlpha(para.mu_mean_alpha, para.mu_sd_alpha);
std::normal_distribution<> normBeta(para.mu_mean_beta, para.mu_sd_beta);
std::normal_distribution<> normD0(para.mu_mean_D0, para.mu_sd_D0);
std::normal_distribution<> normAlphaD(para.mu_mean_alphaD, para.mu_sd_alphaD);
std::normal_distribution<> normBetaD(para.mu_mean_betaD, para.mu_sd_betaD);
std::normal_distribution<> normGamma(para.mu_mean_gamma, para.mu_sd_gamma);


//Function declarations
const string Int2Str(const int x);
const string Float2Str(const double x);
void RunModel(void);
void Initialisation(int xmax, int ymax);
void env_stoch(std::normal_distribution<>);
void local_extinction(void);
void resource_competition(int x, int y, bool adult);
void inheritance(Individual *pup, Individual mom, Individual dad);
void reproduction(void);
void mutations(void);
void dispersal(void);
void survival(void);
void outPar(void);
void outInd_header(void);
void outPop_header(void);

