#ifndef MONOMER_H
#define MONOMER_H

#include "SimulationParameters.h"
#include "Vector.h"
using namespace std;

class Monomer
{
	public:
		Monomer();			
		Monomer(double x, double y, double z);			
	    ~Monomer();
	        
	    void setPosition(double x, double y, double z);	
	        
        double getX(){return Position.getX();}
        double getY(){return Position.getY();}
        double getZ(){return Position.getZ();}
        
        double operator - (Monomer &);				//Returns Squared Distance between arguments//	
        double operator + (Monomer &);				//Returns Lennard Jonnes Potential, returns "0" if distance > CUTOFF//
        double operator * (Monomer &);				//Returns FENE Potential//
	        	
	private:
		Vector Position;
};		

#endif