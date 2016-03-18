#ifndef HELPERFUNCTIONS_H
#define HELPERFUNCTIONS_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>
#include <fstream>
#include "HistogramExtended.h"
#include "Aggregate.h"
using namespace std;

void initializeTemperatures(int commSize,int myRank, double &localTemperature, double temperatureArray[]);

int replicaExchangeUpdate(int commSize, int myRank, MPI_Comm comm, Aggregate &myAggregate, double temperatureArray[]);

double optimizeDisplacementValues(int myRank, Aggregate &myAggregate, double myTemperature);

void printEnergyHistograms(int myRank, ExtendedHistogram &myHistogram, double myTemperature, int fileIndex);


#endif