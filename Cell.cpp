#include "Cell.h"

// Class constructor
Cell::Cell(int xx, int yy)
{
	x = xx;
	y = yy;
	Res = 0.0; 
	N = 0; 
	Nfem = 0;
	Nmale = 0;
	Noffs = 0;
	Foffs = 0;
	Moffs = 0;
	emigrants = 0;
	prop_emigrants = 0.0;
	eps = 0.0;
	sum_BC = 0.0;

}
//----------------------------------------------------------
Cell::~Cell()
{
}
//----------------------------------------------------------
void Cell::DeleteAdults(void){
	std::vector<Individual>::iterator iter;

	//Delete adult females
	iter = females.begin();
	while (iter != females.end()){
		iter->deleteInd();
		iter++;
	}
	//Delete adult males
	iter = males.begin();
	while (iter != males.end()){
		iter->deleteInd();
		iter++;
	}
	females.clear();
	males.clear();

	N = 0;
	Nfem = 0;
	Nmale = 0;
}