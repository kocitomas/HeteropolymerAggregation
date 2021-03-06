//Aggregate Class to be used in Wang-Landau, Metropolis, MUCA, and Parallel-Tempering simulations//
#ifndef AGGREGATE_H
#define AGGREGATE_H

#include <cstring>
#include <fstream>
#include <string>
#include <set>
#include <vector>
#include "Polymer.h"
using namespace std;

class Aggregate
{
	public:
		
		Polymer 	**PolymerArray;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                               // DESIGNATED CONSTRUCTORS //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Constructor for HOMOPOLYMER simulations
		Aggregate(int aggregateSize, int polymerChainLength, double monomerTypeArrayIntra[], double monomerTypeArrayInter[], bool bondedLJ);
		
		~Aggregate();

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                               // GETTERS/SETTERS //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		double getEnergyOfPolymer(int whichPolymer){return PolymerPotentialEnergies[whichPolymer];}
		double getInteractionEnergy(){return AggregateInteractionEnergy;}
		double getTotalEnergy(){return AggregateEnergy;}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                               // DISPLACEMENT UPDATE //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		RESULT updatePositionWithSphericalBoundaries(int whichPolymer,int whichMonomer,double dx, double dy, double dz, double constraintRadiusSquared);
	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                               // ENERGY CALCULATION//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Use this method after each sucessful replica exchange 
		// Recalculates the potential energies of individual polymers as well as the interaction energy of the aggregate
		void recalculateAggregateEnergy();

		// To be called by the designated initializer in order to perform the initial setup of the interaction energy matrix 
		void calculateInteractionEnergy();

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                               // ENERGY CHANGE CALCULATION//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Use with Single Displacement Metropolis Update //
		double findTotalEnergyDifferenceDueToUpdate(int whichPolymer, int whichMonomer);
		


		double findInteractionEnergyDifferenceDueToUpdate(int whichPolymer, int whichMonomer);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                               // ENERGY UPDATE//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Use whenever a displacement updated is accepted 
		void updatePolymerEnergy(int whichPolymer,int whichMonomer);

		void updateTotalEnergy(int whichPolymer, int whichMonomer);
		void updateInteractionEnergy(int whichPolymer,int whichMonomer);
		
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                               // MONTE CARLO UPDATE //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Single Displacement Metropolis Update //
		int performSingleMetropolisUpdate(double constraintRadiusSquared, double updateDisplacementMagnitude, double canonicalTemperature, double maxAllowedEnergy);

		// Global Displacement Metropolis Update //
		int performGlobalMetropolisUpdate(double constraintRadiusSquared, double updateDisplacementMagnitude, double canonicalTemperature, double maxAllowedEnergy);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                               // AGGREGATE STRUCTURAL QUANTITIES //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void calculateCenterOfMass();

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                               // STRUCTURE OUTPUT //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		void printPositionsToFileInMolFormat(ofstream& outputFile);

	private:

		int 		AggregateSize;
		int 		PolymerChainLength;
		double 		AggregateInteractionEnergy;
		double 		AggregateEnergy;

		double 		*InteractionEnergyMatrix;
		double 		*TempInteractionEnergyArray;
		double 		*PolymerPotentialEnergies;

		double 		CenterOfMassX, CenterOfMassY, CenterOfMassZ;


		// Helper variables //
		double InteractionEnergyDifference;
		double PolymerEnergyDifference;

};
#endif