#include "Histogram.h"
using namespace std;


//#################################### CONSTRUCTOR #####################################//
Histogram::Histogram(int numberOfBins, double minValue, double maxValue, double initialValue):max(maxValue),min(minValue),numberOfBins(numberOfBins)
{
	// Initialization of Member Variables //
	binWidth		= (max - min)/numberOfBins;
	pHistogramArray	= new double[numberOfBins];

	//Initialization of Histogram_Array//
	for (int i = 0; i < numberOfBins; i++) pHistogramArray[i] = initialValue;
}

//#################################### DESTRUCTOR #####################################//	
Histogram::~Histogram()
{
	delete [] pHistogramArray;
}

//################################## RESET HISTOGRAM ##################################//
void Histogram::resetHistogramTo(double resetValue)
{
	for (int i = 0; i < numberOfBins; i++) pHistogramArray[i] = resetValue;
}

//################################# UPDATE HISTOGRAM #################################//				
void Histogram::updateAtValueByIncrement(double value, double increment)
{
	int index  = (int)((value - min)/binWidth);
	
	if( value < max && value > min)
	{
		pHistogramArray[index] += increment;
	}
	else
	{
		cout << value << " Value out of range!\n";
	}
}

void Histogram::updateAtIndexByIncrement(int index, double increment)
{
	if(index >= numberOfBins)
	{
		cout << "Index out of range!\n";
	}
	else
	{
		pHistogramArray[index] += increment;
	}
}

void Histogram::setValueAtIndexTo(int index, double value)
{
	pHistogramArray[index] = value;
}

//################################# GET VALUE ########################################//		

double Histogram::getValueAtValue(double value)
{	
	int index		= (int)((value - min)/binWidth);
	return pHistogramArray[index];
}

double Histogram::getValueAtIndex(int index)
{
	return pHistogramArray[index];
}	

		
		
		
		
		
		
		
		
