//Generic Histogram Base Class - Derived classes can extend functionality//
#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <iostream>
#include <cmath>
using namespace std;


class Histogram
{
	public:
		//######################## DESIGNATED CONSTRUCTOR/DESTRUCTOR ##########################//
		//Initializes the histogram with "numberOfBins" bins in the range (min_value,max_value)//																			
		Histogram(int numberOfBins, double minValue, double maxValue, double initialValue);				
		~Histogram();
		
		//############################## RESET ######################################//
		void   resetHistogramTo(double resetValue);		

		//############################# SETTERS #####################################//					        
		void   updateAtValueByIncrement(double value, double increment);					
		void   updateAtIndexByIncrement(int index, double increment);
		void   setValueAtIndexTo(int index, double value);

		//############################# GETTERS #####################################//
		double getValueAtValue(double value);												
		double getValueAtIndex(int index);										
		
		double getBinWidth()const{return binWidth;}
		double getMax()const{return max;}
		double getMin()const{return min;}
		int getNumberOfBins()const{return numberOfBins;}

			
		
	protected:
		double *pHistogramArray;												
		double max, min;
		double binWidth;
		int    numberOfBins;
};		

#endif