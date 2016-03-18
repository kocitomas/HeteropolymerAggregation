#ifndef POLYMER_H 
#define POLYMER_H

#include "Monomer.h"
#include "RandomNumberGenerator.h"
using namespace std;

enum RESULT{ACCEPTED, REJECTED};

class Polymer
{
	public:
		Monomer *MonomerArray;
		double	*MonomerTypeArrayIntra;
		double 	*MonomerTypeArrayInter;
	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                               // CONSTRUCTOR/DESTRUCTOR //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
		
		//Initiates Polymer with N monomers in a straight line with spacing given by R_0
		//The position of the first monomer in the chain is given by (xInit,yInit,zInit)
		
		Polymer(int polymerLength, double xInit, double yInit, double zInit, double monomerTypeArrayIntra[], double monomerTypeArrayInter[], double bondedLJ);
		~Polymer();

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                               		// POLYMER PARAMETERS //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

		int    getPolymerLength()const{return PolymerLength;}
	    double getPolymerEnergy()const{return PolymerEnergy;}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                               		// POSITIONAL UPDATE //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////   
	    
	    // Updates the position a given monomer (indicated by "index") by dx,dy,dz. 
	    // Returns ACCEPTED if the update is legal according to provided spherical constraints  
	    // Note: Enter the squared radius of the spherical constraint volume    
	    
	    RESULT updatePositionWithSphericalBoundaries(int monomerIndex, double dx, double dy, double dz, double constraintRadiusSquared); 

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                               // POTENTIAL ENERGY CALCULATIONS //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

	    // Calculates and returns the current potential energy of the polymer
		// Use whenever the monomer coordinates are updated without modifying the inter-monomer energies
	    void recalculatePolymerEnergy();

	    // Sums the InterMonomerPotentialEnergyMatrix and returns the current total potential energy of the polymer  																
	    double sumPolymerPotentialEnergy() const;

	    // Returns the energy difference between the updated and the old state provided the index of the displaced monomer
	    double getEnergyDifferenceDueToUpdate(int monomerIndex);	

	    // Updates the potential energy array with the data from the temporary array. Use after a proposed update is accepted
	    void   updateInterMonomerEnergyMatrix(int monomerIndex);

	    	  
	private:
		bool	BondedLJ;
		int 	PolymerLength;
		double  PolymerEnergy;
	  	double  *InterMonomerPotentialEnergyMatrix;
	  	double  *MonomerTemporaryPotentialEnergy;
	  	double 	*HeteroPolymerModifiers;
};

#endif