#ifndef SIMULATION_PARAMETERS_H
#define SIMULATION_PARAMETERS_H

#include <iostream>
#include <cmath>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
using namespace std;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                               		// SIMULATION PARAMETERS //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Number of independent processes //
const int NUMBER_OF_REPLICAS        = 92;

// Model parameters //
const double R_0					= 0.7;
const double R						= 0.3;
const double K						= 40.0;
const double FENE					= ((K/2.0)*(R*R));

const double SIGMA					= (R_0/pow(2.0,(1.0/6.0)));
const double SIGMA_2 				= SIGMA*SIGMA;
const double SIGMA_6     	 		= 4.0*(SIGMA_2*SIGMA_2*SIGMA_2);
const double SIGMA_12 				= SIGMA_6*(SIGMA_2*SIGMA_2*SIGMA_2);
const double CUTOFF_POTENTIAL		= 0.016316891136;
const double CUTOFF_SQUARED			= (6.25*SIGMA_2);
const double FENE_TEST_MAX			= ((R_0 + R) * (R_0 + R))-0.001;
const double FENE_TEST_MIN			= ((R_0 - R) * (R_0 - R))+0.001; 

// Aggregate parameters //
const int AGGREGATE_SIZE			= 3;
const int POLYMER_LENGTH     		= 5;
const bool BONDED_LJ				= true;
const double MONOMER_MOD_INTRA		= 1.0;
const double MONOMER_MOD_INTER		= 1.0;
const double PARTICLE_DENSITY		= 0.001;
const double SPHERICAL_CONSTRAINT	= pow(((AGGREGATE_SIZE*POLYMER_LENGTH)/((4.0/3.0)*3.14*PARTICLE_DENSITY)),(2.0/3.0));
const double SINGLE_DISPLACEMENT   	= (R_0/5.0);
const double GLOBAL_DISPLACEMENT    = pow(SPHERICAL_CONSTRAINT,0.5)/20.0;

// Number of iterations and frequency of replica exchange steps //
const int TRAINING_PHASE			= 1000;
const int EQUILIBRATION				= 10000;

// Total Number of Sweeps = DATA_POINTS * REPLICA_EXCHANGES * GLOBAL_UPDATE_FREQ * METROPOLIS_SWEEPS
const int DATA_POINTS				= 10;
const int REPLICA_EXCHANGES		  	= 10000;
const int GLOBAL_UPDATE_FREQ 		= 10;
const int METROPOLIS_SWEEPS     	= (50*POLYMER_LENGTH*AGGREGATE_SIZE);

// Histogram parameters //
const double MIN_ENERGY				= -51;
const double MAX_ENERGY				= 20;
const int NUMBER_OF_BINS 			=5000;

// File output parameters //
const int PRINT_EVERY					= 1000;

#endif
