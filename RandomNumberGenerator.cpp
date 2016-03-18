#include <iostream>
#include <cmath>
#include "RandomNumberGenerator.h"

using namespace std;


void genrandInit(unsigned long long seed)
{
	init_genrand64(seed);
}

double genrandBasic()
{
	return genrand64_real1();
}

int genrandIndex(const int polymer_length)
{
	return (int)(genrand64_real1() * polymer_length);
}

void genrandDisplacementUpdate(const double adjust_factor, double &d_x, double &d_y, double &d_z)
{
	d_x		= (genrand64_real1() - 0.5) * adjust_factor;
	d_y		= (genrand64_real1() - 0.5) * adjust_factor;
	d_z		= (genrand64_real1() - 0.5) * adjust_factor;
}	
			
