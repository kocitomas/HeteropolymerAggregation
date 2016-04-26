#include "Aggregate.h"
#include "HelperFunctions.h"
#include <fstream>
using namespace std;


int main()
{
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                           // MPI Initialization //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	int myRank, commSize;
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &commSize);  

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                      // Temperature Initialization //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	double temperatureArray[NUMBER_OF_REPLICAS];
	double localTemperature;
	initializeTemperatures(commSize,myRank,localTemperature,temperatureArray);

	double aggregateEnergy;
	double aggregatePotentialEnergy;
	double aggregateKineticEnergy;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                           // Output Files //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//char configurationOutputFileName[128];
	//sprintf(configurationOutputFileName,"configAgg_%d_Len_%d_Temp_%f.mol2",AGGREGATE_SIZE,POLYMER_LENGTH,localTemperature);
	//ofstream configurationOutputFile(configurationOutputFileName);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                      // Monomer-Monomer Interaction Strength //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	double monomerTypeArrayIntra[POLYMER_LENGTH];
	monomerTypeArrayIntra[0] = MONOMER_MOD_INTRA;
	for(int i = 1; i < POLYMER_LENGTH; i ++) monomerTypeArrayIntra[i] = 1.0;

	double monomerTypeArrayInter[POLYMER_LENGTH];
	monomerTypeArrayInter[0] = MONOMER_MOD_INTER;
	for(int i = 1; i < POLYMER_LENGTH; i ++) monomerTypeArrayInter[i] = 1.0;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                      // Histogram Initialization //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	ExtendedHistogram *totalEnergyHistogram 		= new ExtendedHistogram(NUMBER_OF_BINS,MIN_ENERGY,MAX_ENERGY,0.0);
	ExtendedHistogram *potentialEnergyHistogram 	= new ExtendedHistogram(NUMBER_OF_BINS,MIN_ENERGY,MAX_ENERGY,0.0);
	ExtendedHistogram *kineticEnergyHistogram 		= new ExtendedHistogram(NUMBER_OF_BINS,MIN_ENERGY,MAX_KINETIC_ENERGY,0.0);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                        // Equilibration and optimization //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	Aggregate *myAggregate = new Aggregate(AGGREGATE_SIZE,POLYMER_LENGTH,monomerTypeArrayIntra,monomerTypeArrayInter,BONDED_LJ);

	// Equilibration-1
	for(int i = 0; i < EQUILIBRATION; i++)
	{
		for(int j = 0; j < METROPOLIS_SWEEPS; j++)
		{
			if(myRank != 0)
			{
				myAggregate -> performSingleMetropolisUpdate(SPHERICAL_CONSTRAINT,SINGLE_DISPLACEMENT,localTemperature,MAX_ENERGY);	 
				if(MOMENTUM_UPDATE)
				{
					myAggregate -> performSingleMomentumMetropolisUpdate(MOMENTUM_DISPLACEMENT,localTemperature,MAX_ENERGY);
				}
			} 
		} 
		replicaExchangeUpdate(commSize,myRank,MPI_COMM_WORLD,*myAggregate,temperatureArray);
	}

	// Displacement optimization
	double trainedDisplacement = 1.0;

	if(myRank != 0)
	{
		trainedDisplacement = optimizeDisplacementValues(myRank, *myAggregate, localTemperature);
		cout << "Process " << myRank << " at temperature " << localTemperature << " displacement optimized at "<<trainedDisplacement <<"\n";
	}
	
	// Equilibration-2
	for(int i = 0; i < EQUILIBRATION; i++)
	{
		for(int j = 0; j < METROPOLIS_SWEEPS; j++)
		{
			if(myRank != 0)
			{
				myAggregate -> performSingleMetropolisUpdate(SPHERICAL_CONSTRAINT,trainedDisplacement,localTemperature,MAX_ENERGY);	 
				if(MOMENTUM_UPDATE)
				{
					myAggregate -> performSingleMomentumMetropolisUpdate(MOMENTUM_DISPLACEMENT,localTemperature,MAX_ENERGY);
				}
			} 
		} 
		replicaExchangeUpdate(commSize,myRank,MPI_COMM_WORLD,*myAggregate,temperatureArray);
	}
	MPI_Barrier(MPI_COMM_WORLD);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                           // Acceptance Ratios //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	double singleDisplacementAcceptanceRatio 	= 0.0;
	double globalDisplacementAcceptanceRatio 	= 0.0;
	double replicaExchangeAcceptanceRatio 		= 0.0;
	double momentumUpdateAcceptanceRatio		= 0.0;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                           // Data Collection //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	
	for (int dataPoint = 0; dataPoint < DATA_POINTS; dataPoint ++)
	{
		// Temporary Histogram //
		ExtendedHistogram *temporaryEnergyHistogram 	= new ExtendedHistogram(NUMBER_OF_BINS,MIN_ENERGY,MAX_ENERGY,0.0);

		for(int i = 0; i < REPLICA_EXCHANGES; i++)
		{
			for(int j = 0; j < METROPOLIS_SWEEPS; j++)
			{
				for(int k = 0; k < GLOBAL_UPDATE_FREQ; k++)
				{
					if(myRank != 0)
					{
						singleDisplacementAcceptanceRatio 	+=   myAggregate 	-> performSingleMetropolisUpdate(SPHERICAL_CONSTRAINT,trainedDisplacement,localTemperature,MAX_ENERGY);
							
						if(MOMENTUM_UPDATE)
						{
							momentumUpdateAcceptanceRatio		+= myAggregate 		-> performSingleMomentumMetropolisUpdate(MOMENTUM_DISPLACEMENT,localTemperature,MAX_ENERGY);
						}

						aggregateKineticEnergy				 = myAggregate		-> getKineticEnergy();
						kineticEnergyHistogram 									-> updateAtValueByIncrement(aggregateKineticEnergy,1.0);
						
						aggregatePotentialEnergy			 = myAggregate		-> getPotentialEnergy();
						potentialEnergyHistogram								-> updateAtValueByIncrement(aggregatePotentialEnergy,1.0);
						
						aggregateEnergy 					 = myAggregate 		-> getTotalEnergy();
						totalEnergyHistogram									-> updateAtValueByIncrement(aggregateEnergy,1.0);
						temporaryEnergyHistogram    							-> updateAtValueByIncrement(aggregateEnergy,1.0);
						//cout << aggregateEnergy <<"\t" << aggregatePotentialEnergy+aggregateKineticEnergy << "\t" << aggregatePotentialEnergy << "\t" << aggregateKineticEnergy << "\t" << "\n";
					} 
				}
				
				if(GLOBAL_UPDATE)
				{
					if(myRank != 0)
					{
						globalDisplacementAcceptanceRatio 	+=   myAggregate 	-> performGlobalMetropolisUpdate(SPHERICAL_CONSTRAINT,GLOBAL_DISPLACEMENT,localTemperature,MAX_ENERGY);
						aggregateEnergy 					 =   myAggregate 	-> getTotalEnergy();
						totalEnergyHistogram									-> updateAtValueByIncrement(aggregateEnergy,1.0);
						temporaryEnergyHistogram    							-> updateAtValueByIncrement(aggregateEnergy,1.0);
					}
				} 
			}

			// Perform Replica Exchange Update //
			replicaExchangeAcceptanceRatio += replicaExchangeUpdate(commSize,myRank,MPI_COMM_WORLD,*myAggregate,temperatureArray);

			if(i%PRINT_EVERY == 0)
			{
				cout << "Process " << myRank <<" reached iteration "<< i << "\n";
				//myAggregate -> printPositionsToFileInMolFormat(configurationOutputFile);
			}	
		}
		printEnergyHistograms(myRank,*temporaryEnergyHistogram, localTemperature, dataPoint);
		delete temporaryEnergyHistogram;
	}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                                      // File Write & Clean-Up //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	singleDisplacementAcceptanceRatio = singleDisplacementAcceptanceRatio/(REPLICA_EXCHANGES*METROPOLIS_SWEEPS*GLOBAL_UPDATE_FREQ*DATA_POINTS);
	globalDisplacementAcceptanceRatio = globalDisplacementAcceptanceRatio/(REPLICA_EXCHANGES*METROPOLIS_SWEEPS*DATA_POINTS);
	replicaExchangeAcceptanceRatio    = replicaExchangeAcceptanceRatio/(REPLICA_EXCHANGES*DATA_POINTS);
	momentumUpdateAcceptanceRatio 	  = momentumUpdateAcceptanceRatio/(REPLICA_EXCHANGES*METROPOLIS_SWEEPS*GLOBAL_UPDATE_FREQ*DATA_POINTS);
	cout << localTemperature << "\t" << singleDisplacementAcceptanceRatio << "\t" << globalDisplacementAcceptanceRatio << "\t" << replicaExchangeAcceptanceRatio << "\t" << momentumUpdateAcceptanceRatio <<"\n";
	     
	printEnergyHistograms(myRank,*totalEnergyHistogram, localTemperature, DATA_POINTS);
	printEnergyHistograms(myRank,*potentialEnergyHistogram, localTemperature, DATA_POINTS+1);
	printEnergyHistograms(myRank,*kineticEnergyHistogram, localTemperature, DATA_POINTS+2);

	printCanonicalQuantitiesfromHistogram(myRank, *totalEnergyHistogram, localTemperature, 1);
	printCanonicalQuantitiesfromHistogram(myRank, *potentialEnergyHistogram, localTemperature, 2);
	printCanonicalQuantitiesfromHistogram(myRank, *kineticEnergyHistogram, localTemperature, 3);
	//configurationOutputFile.close();
	MPI_Finalize(); 
	
	return 0;
}
