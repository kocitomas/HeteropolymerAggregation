#include "Polymer.h"
using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                               	// CONSTRUCTOR/DESTRUCTOR //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Polymer::Polymer(int polymerLength, double xInit, double yInit, double zInit, double monomerTypeArrayIntra[], double monomerTypeArrayInter[], bool bondedLJ):PolymerLength(polymerLength)
{
	// Check for bonded LJ
	BondedLJ = bondedLJ;
	// Create an initialize an array to hold monomer types (inter/intra)
	MonomerTypeArrayIntra = new double[PolymerLength];
	for(int i = 0; i < PolymerLength; i++) MonomerTypeArrayIntra[i] = monomerTypeArrayIntra[i];

	MonomerTypeArrayInter = new double[PolymerLength];
	for(int i = 0; i < PolymerLength; i++) MonomerTypeArrayInter[i] = monomerTypeArrayInter[i];

	// Create an array of monomers and initialize their coordinates to (0,0,0)
	MonomerArray = new Monomer[PolymerLength];								
	
	// Set the positions of the individual monomers to a straight chain with bond lengths of R_0
	for (int i = 0; i < PolymerLength; i ++) MonomerArray[i].setPosition((R_0 * i) + xInit, yInit, zInit); 	
	
	//Create and initialize to "0" matrix which will hold inter-monomer potentials
	InterMonomerPotentialEnergyMatrix = new double[PolymerLength * PolymerLength];

	for (int i = 0; i < PolymerLength; i++)
	{
		for (int j = 0; j < PolymerLength; j++)
		{
			InterMonomerPotentialEnergyMatrix[(i*PolymerLength) + j] = 0.0;
		}
	}
	
	//Update the inter-monomer potential energy matrix to current values
	for (int i = 0; i < PolymerLength - 1; i++)										
	{
		for (int j = i + 1; j < PolymerLength; j++)
		{
			// Check for Lennard-Jonnes potential only since the bond lengths are all set to R_0
			InterMonomerPotentialEnergyMatrix[(i*PolymerLength) + j] = (MonomerArray[i] + MonomerArray[j])*(MonomerTypeArrayIntra[i]*MonomerTypeArrayIntra[j]);
			InterMonomerPotentialEnergyMatrix[(j*PolymerLength) + i] = InterMonomerPotentialEnergyMatrix[(i*PolymerLength) + j];
		}
	}
	
	//Create and initialize the temporary energy array (used with metropolis updates)
	MonomerTemporaryPotentialEnergy = new double[PolymerLength];								
	for (int i = 0; i < PolymerLength; i++) MonomerTemporaryPotentialEnergy[i] = 0.0;
	
	//Create and initalize the polymer momentum array
	MomentumArray = new double[3*PolymerLength];
	for (int i = 0; i < 3*PolymerLength; i++) MomentumArray[i] = 0.0;

	//Calculate polymer potential energy
	PolymerPotentialEnergy = sumPolymerPotentialEnergy();
	PolymerKineticEnergy   = 0.0;
	PolymerTotalEnergy	   = PolymerPotentialEnergy + PolymerKineticEnergy;
}

Polymer::~Polymer()
{
	delete [] MonomerArray;
	delete [] InterMonomerPotentialEnergyMatrix;
	delete [] MonomerTemporaryPotentialEnergy;
	delete [] MonomerTypeArrayIntra;
	delete [] MomentumArray;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                               		// POSITIONAL UPDATE //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

RESULT Polymer::updatePositionWithSphericalBoundaries(int monomerIndex, double dx, double dy, double dz, double constraintRadiusSquared)
{
	double xUpdated = MonomerArray[monomerIndex].getX() + dx;
	double yUpdated = MonomerArray[monomerIndex].getY() + dy;
	double zUpdated = MonomerArray[monomerIndex].getZ() + dz;
	
	// Update monomer position
	MonomerArray[monomerIndex].setPosition(xUpdated, yUpdated, zUpdated);
	
	// Check if updated position is within spherical constraint
	double distanceSquaredFromOrigin = (xUpdated*xUpdated + yUpdated*yUpdated + zUpdated*zUpdated);
	if (distanceSquaredFromOrigin >= constraintRadiusSquared)
	{
		return REJECTED;
	}

	return ACCEPTED;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                               		// MOMENTUM UPDATE //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Polymer::updateMomentum(int monomerIndex, double dPx, double dPy, double dPz)
{
	MomentumArray[(3*monomerIndex)] 	+= dPx;
	MomentumArray[(3*monomerIndex)+1] 	+= dPy;
	MomentumArray[(3*monomerIndex)+2] 	+= dPz;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                               // ENERGY CALCULATIONS //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
void Polymer::calculatePolymerTotalEnergy()
{
	calculatePolymerPotentialEnergy();
	calculatePolymerKineticEnergy();

	PolymerTotalEnergy = PolymerPotentialEnergy + PolymerKineticEnergy;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                               // POTENTIAL ENERGY CALCULATIONS //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

void Polymer::calculatePolymerPotentialEnergy()
{
	// Reset energy array to zero
	for (int i = 0; i < PolymerLength; i++)
	{
		for (int j = 0; j < PolymerLength; j++)
		{
			InterMonomerPotentialEnergyMatrix[(i*PolymerLength) + j] = 0.0;
		}
	}
	
	// Update the inter-monomer potential energy matrix to current energy values	
	for (int i = 0; i < PolymerLength - 1; i++)										
	{
		for (int j = i + 1; j < PolymerLength; j++)
		{
			// Bonded Potential
			if(j == i+1)
			{
				// LJ
				if(BondedLJ) 
				{
					InterMonomerPotentialEnergyMatrix[(i*PolymerLength) + j] += (MonomerArray[i] + MonomerArray[j])*(MonomerTypeArrayIntra[i]*MonomerTypeArrayIntra[j]);
				}
				// FENE
				InterMonomerPotentialEnergyMatrix[(i*PolymerLength) + j] += (MonomerArray[i] * MonomerArray[j]);
			} 
			// Non-Bonded Potential
			else 
			{
				InterMonomerPotentialEnergyMatrix[(i*PolymerLength) + j] += (MonomerArray[i] + MonomerArray[j])*(MonomerTypeArrayIntra[i]*MonomerTypeArrayIntra[j]);
			}

			InterMonomerPotentialEnergyMatrix[(j*PolymerLength) + i]  = InterMonomerPotentialEnergyMatrix[(i*PolymerLength) + j];
		}
	}

	// Update the PolymerPotentialEnergy member variable
	PolymerPotentialEnergy = sumPolymerPotentialEnergy();
}
	
double Polymer::sumPolymerPotentialEnergy() const
{	
	double polymerPotentialEnergy	= 0.0;
	
	// Sum inter-monomer potential energy matrix to get the total polymer energy
	for (int i = 0; i < (PolymerLength - 1); i++)
	{
		for(int j = (i+1); j < PolymerLength; j++)
		{
			polymerPotentialEnergy += InterMonomerPotentialEnergyMatrix[(i*PolymerLength)+j];
		}
	}
	
	return polymerPotentialEnergy;
}

double Polymer::getEnergyDifferenceDueToPositionalUpdate(int monomerIndex)
{
	double energyDifference = 0.0;
	
	for (int i = 0; i < PolymerLength; i ++)
	{
		if (i == monomerIndex)
		{
			MonomerTemporaryPotentialEnergy[i] = 0.0;
		}
		else if ((i == (monomerIndex + 1)) || (i == (monomerIndex - 1)))
		{
			// FENE
			MonomerTemporaryPotentialEnergy[i]	 = (MonomerArray[monomerIndex] * MonomerArray[i]);

			// LJ
			if(BondedLJ)
			{
				MonomerTemporaryPotentialEnergy[i] += ((MonomerArray[monomerIndex] + MonomerArray[i])*(MonomerTypeArrayIntra[monomerIndex]*MonomerTypeArrayIntra[i]));
 			}
			 
			energyDifference +=  MonomerTemporaryPotentialEnergy[i] - InterMonomerPotentialEnergyMatrix[(monomerIndex*PolymerLength) + i];
		}
		else
		{
			MonomerTemporaryPotentialEnergy[i] 	 = (MonomerArray[monomerIndex] + MonomerArray[i]);
			energyDifference					+=  MonomerTemporaryPotentialEnergy[i] - InterMonomerPotentialEnergyMatrix[(monomerIndex*PolymerLength) + i];	
		}
	}
	return energyDifference;
}						
		
void Polymer::updateInterMonomerEnergyMatrix(int monomerIndex)
{	
	for (int i = 0; i < PolymerLength; i++)
	{
		InterMonomerPotentialEnergyMatrix[(monomerIndex*PolymerLength) + i] = MonomerTemporaryPotentialEnergy[i];
		InterMonomerPotentialEnergyMatrix[(i*PolymerLength) + monomerIndex] = MonomerTemporaryPotentialEnergy[i];
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                               // KINETIC ENERGY CALCULATIONS //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 		

void Polymer::calculatePolymerKineticEnergy()
{
	double kineticEnergy = 0.0;

	for(int i = 0; i < 3*PolymerLength; i ++)
	{
		kineticEnergy += (MomentumArray[i]*MomentumArray[i])/2.0;
	}

	PolymerKineticEnergy = kineticEnergy;
}

double Polymer::getKineticEnergyOfMonomer(int whichMonomer)
{
	double kineticEnergy = 0.0;

	for(int i = 0; i < 3; i ++)
	{
		kineticEnergy += (MomentumArray[(whichMonomer*3)+i]*MomentumArray[(whichMonomer*3)+i])/2.0;
	}

	return kineticEnergy;
}