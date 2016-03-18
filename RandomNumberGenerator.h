// Functions for generating random numbers in convenient form //
// Don't forget to seed the random number generator if you want a different sequence of random numbers for each simulation //
#ifndef RANDOM_NO_H
#define RANDOM_NO_H

#include <iostream>
#include <cmath>
#include "MersenneTwister.h"

using namespace std;


// Initialize the random number generator with a seed
void genrandInit(unsigned long long seed);

// Generates a random double in the interval (0,1) //
double genrandBasic();

// Generates a random integer on the interval [0,max_index - 1]
int genrandIndex(const int max_index);

// Generates a displacement vector with the length given by magnitude //
// The vector components are to be passed by reference //
void genrandDisplacementUpdate(const double magnitude, double &, double &, double &);


#endif