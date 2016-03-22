//Aggregate Class to be used in Wang-Landau, Metropolis, MUCA, and Parallel-Tempering simulations//

#include "Aggregate.h"
using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                               // DESIGNATED CONSTRUCTOR/DESTRUCTOR //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Aggregate::Aggregate(int aggregateSize, int PolymerChainLength, double monomerTypeArrayIntra[], double monomerTypeArrayInter[], bool bondedLJ):AggregateSize(aggregateSize),PolymerChainLength(PolymerChainLength),InteractionEnergyDifference(0.0),PolymerEnergyDifference(0.0)
{
	// Create array of pointers to polymer objects
	PolymerArray = new Polymer*[AggregateSize];

	// Assign polymer objects to the pointers in the PolymerArray 
	for(int i = 0; i < AggregateSize; i ++) PolymerArray[i] = new Polymer(PolymerChainLength,0.0,i*R_0,0.0,monomerTypeArrayIntra,monomerTypeArrayInter, bondedLJ);
	
	// Store individual polymer energies 
	PolymerPotentialEnergies = new double[AggregateSize];
	for(int i = 0; i < AggregateSize; i++) PolymerPotentialEnergies[i] = PolymerArray[i]-> getPolymerEnergy();

	// Initialize the interaction energy matrix and the temporary interaction energy array 
	InteractionEnergyMatrix 			= new double[(PolymerChainLength*AggregateSize)*(PolymerChainLength*AggregateSize)];
	TempInteractionEnergyArray 			= new double[PolymerChainLength*AggregateSize];	

	for(int i = 0; i < (PolymerChainLength*AggregateSize)*(PolymerChainLength*AggregateSize); i++) InteractionEnergyMatrix[i] = 0.0;
	for(int i = 0; i < (PolymerChainLength*AggregateSize); i ++) TempInteractionEnergyArray[i] = 0.0;

	// Calculate the total interaction energy
	calculateInteractionEnergy();

	// Calculate the total potential energy of the aggregate
	AggregateEnergy = AggregateInteractionEnergy;

	for(int i = 0; i < AggregateSize; i++) AggregateEnergy += PolymerPotentialEnergies[i];	
}

Aggregate::~Aggregate()
{
	for(int i = 0; i < AggregateSize; i++){
		delete PolymerArray[i];
	}
	delete[] PolymerArray;
	delete[] PolymerPotentialEnergies;
	delete[] InteractionEnergyMatrix;
	delete[] TempInteractionEnergyArray;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                               		// POSITIONAL UPDATES //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

RESULT Aggregate::updatePositionWithSphericalBoundaries(int whichPolymer,int whichMonomer,double dx, double dy, double dz, double constraintRadiusSquared)
{
	RESULT result = PolymerArray[whichPolymer]-> updatePositionWithSphericalBoundaries(whichMonomer, dx, dy, dz, constraintRadiusSquared);
	return result;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                               // ENERGY CALCULATION //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Aggregate::recalculateAggregateEnergy()
{
	// Calculate interaction energy
	calculateInteractionEnergy();

	// Calculate polymer energies
	for(int i = 0; i < AggregateSize; i++){
		PolymerArray[i] -> recalculatePolymerEnergy();
		PolymerPotentialEnergies[i] = PolymerArray[i] -> getPolymerEnergy();
	} 

	// Update the Aggregate Energy member variable
	AggregateEnergy = AggregateInteractionEnergy;
	for(int i = 0; i < AggregateSize; i++) AggregateEnergy += PolymerPotentialEnergies[i];	
}

void Aggregate::calculateInteractionEnergy()
{
	int row,column;
	int dimension 				= PolymerChainLength*AggregateSize;
	AggregateInteractionEnergy 	= 0.0;

	// Loop through individual polymers
	for(int i = 0; i < (AggregateSize - 1); i++){

		// Loop through individual monomers
		for(int j = 0; j < PolymerChainLength; j++){

			//Loop through interactions
			for(int k = i+1; k < AggregateSize; k++){
				for(int l = 0; l < PolymerChainLength; l++){
					row = ((i*PolymerChainLength)+j);
					column = ((k)*PolymerChainLength)+l;

					InteractionEnergyMatrix[(row*dimension)+(column)] = ((PolymerArray[i]->MonomerArray[j]) + (PolymerArray[k]->MonomerArray[l]))
																	  * ((PolymerArray[i]->MonomerTypeArrayInter[j]) * (PolymerArray[k]->MonomerTypeArrayInter[l]));
					InteractionEnergyMatrix[(column*dimension) + row] = InteractionEnergyMatrix[(row*dimension)+(column)];
					AggregateInteractionEnergy += InteractionEnergyMatrix[(row*dimension)+(column)];
				}
			}
		}
	}	
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                               // ENERGY CHANGE CALCULATION//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double Aggregate::findTotalEnergyDifferenceDueToUpdate(int whichPolymer, int whichMonomer)
{
	InteractionEnergyDifference = findInteractionEnergyDifferenceDueToUpdate(whichPolymer,whichMonomer);
	PolymerEnergyDifference 	= PolymerArray[whichPolymer] -> getEnergyDifferenceDueToUpdate(whichMonomer);
	return (InteractionEnergyDifference + PolymerEnergyDifference);
}

double Aggregate::findInteractionEnergyDifferenceDueToUpdate(int whichPolymer, int whichMonomer)
{
	double energyDifference 		= 0.0;
	int row							= (whichPolymer*PolymerChainLength) + whichMonomer;

	for(int i = 0; i < AggregateSize; i++){
		for(int j = 0; j < PolymerChainLength; j++){
			if(i == whichPolymer){
				TempInteractionEnergyArray[(i*PolymerChainLength) + j] = 0.0;
			}
			else{
				TempInteractionEnergyArray[(i*PolymerChainLength) + j] = ((PolymerArray[whichPolymer] -> MonomerArray[whichMonomer])+(PolymerArray[i] -> MonomerArray[j]))
																		 * (PolymerArray[whichPolymer]->MonomerTypeArrayInter[whichMonomer] * PolymerArray[i]->MonomerTypeArrayInter[j]);
				energyDifference += (TempInteractionEnergyArray[(i*PolymerChainLength) + j] - InteractionEnergyMatrix[(row*PolymerChainLength*AggregateSize)+(i*PolymerChainLength)+j]);
			}
		}
	}
	return energyDifference;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                               // ENERGY UPDATE//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Aggregate::updateTotalEnergy(int whichPolymer, int whichMonomer)
{
	updateInteractionEnergy(whichPolymer,whichMonomer);
	updatePolymerEnergy(whichPolymer,whichMonomer);
	AggregateEnergy = AggregateInteractionEnergy;

	for(int i = 0; i < AggregateSize; i++) AggregateEnergy += PolymerPotentialEnergies[i];	
}

void Aggregate::updateInteractionEnergy(int whichPolymer,int whichMonomer)
{
	int row 		= (whichPolymer*PolymerChainLength) + whichMonomer;
	int dimension 	= (PolymerChainLength*AggregateSize);

	for(int i = 0; i < dimension; i++)
	{
		InteractionEnergyMatrix[(row*dimension) + i] = TempInteractionEnergyArray[i];
		InteractionEnergyMatrix[(i*dimension) + row] = TempInteractionEnergyArray[i];
	}
	AggregateInteractionEnergy += InteractionEnergyDifference;
}

void Aggregate::updatePolymerEnergy(int whichPolymer,int whichMonomer)
{
	PolymerArray[whichPolymer] -> updateInterMonomerEnergyMatrix(whichMonomer);
	PolymerPotentialEnergies[whichPolymer] += PolymerEnergyDifference;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                               // METROPOLIS SINGLE DISPLACEMENT UPDATE //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int Aggregate::performSingleMetropolisUpdate(double constraintRadiusSquared, double updateDisplacementMagnitude, double canonicalTemperature, double maxAllowedEnergy)
{
	double dx,dy,dz;

	// Choose randomly both the polymer and the monomer to be updated
	int whichPolymer = genrandIndex(AggregateSize);
	int whichMonomer = genrandIndex(PolymerChainLength);

	// Generate a normalized triplet of displacement coordinates
	genrandDisplacementUpdate(updateDisplacementMagnitude,dx,dy,dz);

	// Perform a displacement update within fixed spherical boundaries
	RESULT updateSuccess = updatePositionWithSphericalBoundaries(whichPolymer,whichMonomer, dx,dy,dz, constraintRadiusSquared);

	// Revert to previous configuration and exit the function if the update was not accepted
	if (updateSuccess == REJECTED){
		updateSuccess = updatePositionWithSphericalBoundaries(whichPolymer,whichMonomer, -dx,-dy,-dz, constraintRadiusSquared);
		if(updateSuccess == REJECTED) cout << "fatal error during position update \n";
		return 0;
	}

	// Calculate the energy difference due to the update
	double energyDifference = findTotalEnergyDifferenceDueToUpdate(whichPolymer,whichMonomer);

	// Check for nan!!!
	if(energyDifference != energyDifference)
	{
		updateSuccess = updatePositionWithSphericalBoundaries(whichPolymer,whichMonomer, -dx,-dy,-dz, constraintRadiusSquared);
		if(updateSuccess == REJECTED) cout << "fatal error during position update \n";
		return 0;
	}

	// Accept the update if the energy of the new state is less or equal to that of the old state
	if(energyDifference <= 0.0){
		updateTotalEnergy(whichPolymer,whichMonomer);
		return 1;
	}

	// Check if the total energy is withing the defined bounds or reject the update
	double tempAggregateEnergy = AggregateEnergy + energyDifference;

	if(tempAggregateEnergy >= maxAllowedEnergy){
		updateSuccess = updatePositionWithSphericalBoundaries(whichPolymer,whichMonomer, -dx,-dy,-dz, constraintRadiusSquared);
		if(updateSuccess == REJECTED) cout << "fatal error during position update \n";
		return 0;
	}

	// Calculate the exponential term to be used in the metropolis criterion and generate a random number in the interval (0,1)
	double expFactor = exp(-(energyDifference/canonicalTemperature));
	double randFactor = genrandBasic();

	// Reject the update if the exponential term is less or equal to the randomly generated number in the interval (0,1)
	if(randFactor > expFactor){
		updateSuccess = updatePositionWithSphericalBoundaries(whichPolymer,whichMonomer, -dx,-dy,-dz, constraintRadiusSquared);
		if(updateSuccess == REJECTED) cout << "fatal error during position update \n";
		return 0;
	}
	else{
		updateTotalEnergy(whichPolymer,whichMonomer);
		return 1;
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                               // METROPOLIS GLOBAL DISPLACEMENT UPDATE //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int Aggregate::performGlobalMetropolisUpdate(double constraintRadiusSquared, double updateDisplacementMagnitude, double canonicalTemperature, double maxAllowedEnergy)
{
	RESULT updateSuccess;
	double dx,dy,dz;

	// Choose randomly the polymer to be updated
	int whichPolymer = genrandIndex(AggregateSize);

	// Generate a normalized triplet of displacement coordinates
	genrandDisplacementUpdate(updateDisplacementMagnitude,dx,dy,dz);

	// Holds the energy change due to the global update //
	double globalEnergyDifference = 0.0;

	// Perform a displacement update within fixed spherical boundaries for the whole polymer
	for (int i = 0; i < PolymerChainLength; i ++)
	{
		updateSuccess = updatePositionWithSphericalBoundaries(whichPolymer, i, dx,dy,dz, constraintRadiusSquared);
		
		// Revert to previous configuration and exit the function if the updated position is out of bounds
		if (updateSuccess == REJECTED){
			for (int j = i; j >= 0; j--)
			{
				updateSuccess = updatePositionWithSphericalBoundaries(whichPolymer, j, -dx,-dy,-dz, constraintRadiusSquared);
				if(updateSuccess == REJECTED) cout << "fatal error during position update \n";
			}
			return 0;
		}

		// Calculate the energy difference due to the update
		globalEnergyDifference += findInteractionEnergyDifferenceDueToUpdate(whichPolymer,i);

		// Check for nan!!!
		if(globalEnergyDifference != globalEnergyDifference)
		{
			for (int j = i; j >= 0; j--)
			{
				updateSuccess = updatePositionWithSphericalBoundaries(whichPolymer, j, -dx,-dy,-dz, constraintRadiusSquared);
				if(updateSuccess == REJECTED) cout << "fatal error during position update \n";
			}
			return 0;
		}
	}
	
	// Accept the update if the energy of the new state is less or equal to that of the old state
	if(globalEnergyDifference <= 0.0){
		recalculateAggregateEnergy();
		return 1;
	}

	if((AggregateEnergy + globalEnergyDifference) >= maxAllowedEnergy){
		for (int j = PolymerChainLength-1; j >= 0; j--)
			{
				updateSuccess = updatePositionWithSphericalBoundaries(whichPolymer, j, -dx,-dy,-dz, constraintRadiusSquared);
				if(updateSuccess == REJECTED) cout << "fatal error during position update \n";
			}
		return 0;
	}

	// Calculate the exponential term to be used in the metropolis criterion and generate a random number in the interval (0,1)
	double expFactor = exp(-(globalEnergyDifference/canonicalTemperature));
	double randFactor = genrandBasic();

	// Reject the update if the exponential term is less or equal to the randomly generated number in the interval (0,1)
	if(randFactor > expFactor){
		for (int j = PolymerChainLength-1; j >= 0; j--)
			{
				updateSuccess = updatePositionWithSphericalBoundaries(whichPolymer, j, -dx,-dy,-dz, constraintRadiusSquared);
				if(updateSuccess == REJECTED) cout << "fatal error during position update \n";
			}
		return 0;
	}
	else{
		recalculateAggregateEnergy();
		return 1;
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                               // AGGREGATE STRUCTURAL QUANTITIES //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Aggregate::calculateCenterOfMass()
{
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;

	for(int i = 0; i < AggregateSize; i ++)
	{
		for(int j = 0; j < PolymerChainLength; j++)
		{
			x	+= PolymerArray[i] -> MonomerArray[j].getX();
			y	+= PolymerArray[i] -> MonomerArray[j].getY();
			z	+= PolymerArray[i] -> MonomerArray[j].getZ();
		}
	}
	CenterOfMassX = x/(AggregateSize*PolymerChainLength);
	CenterOfMassY = y/(AggregateSize*PolymerChainLength);
	CenterOfMassZ = z/(AggregateSize*PolymerChainLength);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                               // STRUCTURE OUTPUT //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Aggregate::printPositionsToFileInMolFormat(ofstream& outputFile)
{
	double x,y,z;
	calculateCenterOfMass();

	// 1 - <TRIPOS> MOLECULE
	outputFile<< "@<TRIPOS>MOLECULE" <<"\n";
	outputFile<< "aggregate" <<"\n";
	outputFile<< AggregateSize*PolymerChainLength <<" "<< (AggregateSize*PolymerChainLength - AggregateSize) <<" "<< AggregateSize <<" "<< 0 << "\t" << 0 <<"\n";
	outputFile<< "SMALL" <<"\n";
	outputFile<< "NO_CHARGES" <<"\n";
	
	// 2 - <TRIPOS> ATOM
	outputFile<< "@<TRIPOS>ATOM" <<"\n";
	for(int i = 0; i < AggregateSize; i++)
	{
		for(int j = 0; j < PolymerChainLength; j ++)
		{
			outputFile<< ((i * PolymerChainLength) + j + 1) <<"\t";
			outputFile<< "C" << ((i * PolymerChainLength) + j + 1) <<"\t";

			// 1 - Get the coordinates of the current monomer 
			x	= PolymerArray[i] -> MonomerArray[j].getX();
			y	= PolymerArray[i] -> MonomerArray[j].getY();
			z	= PolymerArray[i] -> MonomerArray[j].getZ();
			// 2 - Perform the change of coordinates to the COM system
			x 	= x - CenterOfMassX;
			y 	= y - CenterOfMassY;
			z 	= z - CenterOfMassZ;

			outputFile<< x <<"\t"<< y << "\t" << z << "\t";

			if((PolymerArray[i]->MonomerTypeArrayIntra[j]) == 1.0)
			{
				outputFile<< "C" <<"\t"<< (i+1) <<"\t"<<"Unk"<<"\t"<< 0 <<"\n";
			}
			else
			{
				outputFile<< "H" <<"\t"<< (i+1) <<"\t"<<"Unk"<<"\t"<< 0 <<"\n";
			}
			
		}
	}

	// 3 - <TRIPOS> BOND
	outputFile<< "@<TRIPOS>BOND" <<"\n";
	for(int i = 0; i < AggregateSize; i++)
	{
		for(int j = 0; j < PolymerChainLength-1; j ++)
		{
			outputFile<< ((i * PolymerChainLength) + j + 1) <<"\t";
			outputFile<< ((i * PolymerChainLength) + j + 1) <<"\t";
			outputFile<< ((i * PolymerChainLength) + j + 2) <<"\t";
			outputFile<< 1 <<"\n";
		}
	}
}
