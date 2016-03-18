#include "Monomer.h"
using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                               		// CONSTRUCTORS //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Monomer::Monomer():Position()
{
}

Monomer::Monomer(double x, double y, double z):Position(x,y,z)
{
}

Monomer::~Monomer()
{
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                               		// SET MONOMER COORDINATES //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Monomer::setPosition(double x, double y, double z)
{
	Position.setX(x);
	Position.setY(y);
	Position.setZ(z);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                               		// LJ POTENTIAL //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Monomer::operator + (Monomer & rhs)
{
	double distanceSquared;
	double potentialEnergy;
	double lj6;
	double lj12;
	
	distanceSquared	= ((Position.getX() - rhs.Position.getX()) * (Position.getX() - rhs.Position.getX()))
				+ ((Position.getY() - rhs.Position.getY()) * (Position.getY() - rhs.Position.getY()))	
				+ ((Position.getZ() - rhs.Position.getZ()) * (Position.getZ() - rhs.Position.getZ()));
	
	if (distanceSquared > CUTOFF_SQUARED) 
	{
		return(0.0);
	}
	else
	{	
		lj6 				= distanceSquared * distanceSquared * distanceSquared;
		lj12				= lj6 * lj6; 	
		potentialEnergy		= (SIGMA_12/lj12) - (SIGMA_6/lj6) + CUTOFF_POTENTIAL;
		return(potentialEnergy);	
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                               		// FENE POTENTIAL //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Monomer::operator * (Monomer & rhs)
{
	double distanceSquared;
	double distance;
	double potentialEnergy;
	
	distanceSquared	= ((Position.getX() - rhs.Position.getX()) * (Position.getX() - rhs.Position.getX()))
				+ ((Position.getY() - rhs.Position.getY()) * (Position.getY() - rhs.Position.getY()))	
				+ ((Position.getZ() - rhs.Position.getZ()) * (Position.getZ() - rhs.Position.getZ()));
	distance		= sqrt(distanceSquared);
	
	potentialEnergy	= -FENE*log(1-(((distance - R_0)/R)*((distance - R_0)/R)));
	
	return(potentialEnergy);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                               		// SQUARED DISTANCE //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Monomer::operator - (Monomer & rhs)
{
	double distanceSquared;
	distanceSquared	= ((Position.getX() - rhs.Position.getX()) * (Position.getX() - rhs.Position.getX()))
				+ ((Position.getY() - rhs.Position.getY()) * (Position.getY() - rhs.Position.getY()))	
				+ ((Position.getZ() - rhs.Position.getZ()) * (Position.getZ() - rhs.Position.getZ()));
	
	return(distanceSquared);
}
