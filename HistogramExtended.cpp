#include "HistogramExtended.h"
using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                                         // Designated Constructors //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//#################################### CONSTRUCTOR ########################################//
ExtendedHistogram::ExtendedHistogram(int numberOfBins, double minValue, double maxValue, double initialValue):Histogram(numberOfBins,minValue,maxValue,initialValue)
{
}

//#################################### NORMALIZE HISTOGRAM ########################################//
void ExtendedHistogram::normalizeHistogram()
{
	double weight	= 0.0;
	
	for(int i = 0; i < numberOfBins; i ++) weight	+= pHistogramArray[i];
	if(weight > 0.0)
	{
		for(int i = 0; i < numberOfBins; i ++) pHistogramArray[i] = pHistogramArray[i]/weight;
	}
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                                         // Statistics //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//#################################### GET AVERAGE ########################################//
double ExtendedHistogram::getAverage()
{
	double average 		= 0.0;
	double weight		= 0.0;
	
	for(int i = 0; i < numberOfBins; i ++)
	{
		average += pHistogramArray[i] * (min + (binWidth * i));
		weight	+= pHistogramArray[i];
	}
	average			= average/weight;
	return average;
}
	
//#################################### GET AVERAGE OF SQUARES ########################################//
double ExtendedHistogram::getAverageOfSquares()
{
	double squaredAverage 		= 0.0;
	double weight				= 0.0;
	
	for(int i = 0; i < numberOfBins; i ++)
	{
		squaredAverage 		+= pHistogramArray[i] * (min + (binWidth * i)) * (min + (binWidth * i));
		weight	 		 	+= pHistogramArray[i];
	}
	
	squaredAverage			= squaredAverage/weight;
	return squaredAverage;
}

//#################################### GET STANDARD DEVIATION ########################################//
double ExtendedHistogram::getStandardDeviation()
{
	double average			= getAverage();
	double squaredAverage	= getAverageOfSquares();	
	double sDeviation		= sqrt(squaredAverage - (average * average));
	return sDeviation;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                                         // Print Methods //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ExtendedHistogram::printHistogram(int index)
{
	FILE* fOut;
	char  fOutName[128];
	sprintf(fOutName,"histogram.%d",index);
	fOut	= fopen(fOutName,"w");

	for(int i = 0; i < numberOfBins; i ++) fprintf(fOut,"%f\t%.12f\n",(min + (binWidth * i)),pHistogramArray[i]);

	fclose(fOut);
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                                    // Special Copy Constructors //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//########################## SINGLE HISTOGRAM REWEIGHTING CONSTRUCTOR #################################//
ExtendedHistogram::ExtendedHistogram(ExtendedHistogram & rhs, int numberOfBins, double minValue, double maxValue, double temperature)
:Histogram(numberOfBins,minValue,maxValue,0.0)
{
	double weight	= 0.0;

	for(int i = 0; i < numberOfBins; i ++) weight += exp(-(min + (binWidth * i))/temperature);
	for(int i = 0; i < numberOfBins; i ++) pHistogramArray[i] = ((rhs.getValueAtIndex(i))*exp(-(min + (binWidth * i))/temperature))/weight;	
}

//############################## DERIVATIVE ################################//		
ExtendedHistogram::ExtendedHistogram(ExtendedHistogram & rhs, int numberOfBins, double minValue, double maxValue, int derivativeOrder)
:Histogram(numberOfBins,minValue,maxValue,0.0)
{
	if(derivativeOrder == 1)
	{
		for(int i = 2; i < numberOfBins-2; i++)
		{
			pHistogramArray[i] = (8*rhs.getValueAtIndex(i+1))+(rhs.getValueAtIndex(i-2))-(rhs.getValueAtIndex(i+2))-(8*rhs.getValueAtIndex(i-1));
			pHistogramArray[i] = pHistogramArray[i]/(12*binWidth);
		}
	}
	else
	{
		for(int i = 2; i < numberOfBins-2; i++)
		{
			pHistogramArray[i] = (16*rhs.getValueAtIndex(i+1))+(16*rhs.getValueAtIndex(i-1))-(rhs.getValueAtIndex(i+2))
			-(rhs.getValueAtIndex(i-2))-(30*rhs.getValueAtIndex(i));
			pHistogramArray[i] = pHistogramArray[i]/(12*binWidth*binWidth);
		}
	}
}

//############################## BEZIERE SMOOTH ################################//
ExtendedHistogram::ExtendedHistogram(ExtendedHistogram & rhs, int leftofset, int rightofset):Histogram(rhs.numberOfBins,rhs.min,rhs.max,0.0)
{
	double logOmega = 0.0;
	double logSigma = 0.0;
	double currentEnergy = 0.0;

	double shiftedmin 			= min + leftofset*binWidth;
	double shiftedMax			= max - rightofset*binWidth;
	int    shiftednumberOfBins 	= numberOfBins - leftofset - rightofset;


	for(int i = 0;i < shiftednumberOfBins; i ++)
	{
		currentEnergy = shiftedmin + i*binWidth;
		logOmega 	  = shiftednumberOfBins*log(shiftedMax - currentEnergy) + log(rhs.getValueAtIndex(leftofset)) 
					    - shiftednumberOfBins*log(shiftedMax - shiftedmin);
		logSigma	  = logOmega;

		for(int j = 0; j < shiftednumberOfBins; j ++)
		{
			logOmega = logOmega + log((shiftednumberOfBins - j)/(j+1.0)) + log((currentEnergy - shiftedmin)/(shiftedMax - currentEnergy)) 
			+ log(rhs.getValueAtIndex(leftofset+j+1)/rhs.getValueAtIndex(leftofset+j));

			if (logOmega > logSigma)
			{
				logSigma = logOmega + log(1+exp(logSigma - logOmega));
			} 
			else
			{
				logSigma = logSigma + log(1+exp(logOmega - logSigma));
			} 		
		}
		pHistogramArray[leftofset+i] = exp(logSigma);
	}
}
