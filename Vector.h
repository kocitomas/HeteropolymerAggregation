//Vector class designed specifically for the use by the Monomer class//
#ifndef VECTOR_H
#define VECTOR_H
#include <iostream>
using namespace std;

class Vector 
{
	public:

		Vector();
		Vector(double x, double y, double z);
	    ~Vector(){}
	       
	    double getX(){return xCord;}
	    double getY(){return yCord;}
	    double getZ(){return zCord;}
	       	
	    void   setX(double x){xCord = x;}
	    void   setY(double y){yCord = y;}
	    void   setZ(double z){zCord = z;}
	       	
	private:
		double xCord;
		double yCord;
		double zCord;
};

#endif