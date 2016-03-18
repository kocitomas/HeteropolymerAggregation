//Enhanced Histogram Base Class//
#ifndef HISTOGRAM_EXTENDED_H
#define HISTOGRAM_EXTENDED_H


#include <iostream>
#include <cmath>
#include "Histogram.h"
#include <stdio.h>
using namespace std;


class ExtendedHistogram: public Histogram
{
	public:
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                                         // Designated Constructors //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		ExtendedHistogram(int numberOfBins, double minValue, double maxValue, double initialValue);
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                                         // Special Copy Constructors //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//single histogram reweighting from density of states//
		ExtendedHistogram(ExtendedHistogram &, int numberOfBins, double minValue, double maxValue, double temperature);

		//derivatives of the histogram via five point method, up to second order//
		ExtendedHistogram(ExtendedHistogram &, int numberOfBins, double minValue, double maxValue, int derivativeOrder);

		//Beziere smoothing//
		ExtendedHistogram(ExtendedHistogram &,int leftofset, int rightofset);
		
		void   normalizeHistogram();
		double getAverage();
		double getAverageOfSquares();
		double getStandardDeviation();


		void   printHistogram(int index);
};	

#endif