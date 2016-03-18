#include "HelperFunctions.h"
using namespace std;

#define MIN(a,b) ((a) < (b) ? (a) : (b))

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                               // TEMPERATURE INITIALIZATION //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void initializeTemperatures(int commSize,int myRank, double &localTemperature, double temperatureArray[])
{
	FILE *temperatures;
	temperatures 			= fopen("temperatures","r");
	temperatureArray[0]    	= 0.0;

	for(int i = 1; i < commSize; i ++) fscanf(temperatures,"%lf\n",&temperatureArray[i]);		

	localTemperature  	= temperatureArray[myRank];
	fclose(temperatures);
} 

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                               // REPLICA EXCHANGE UPDATE //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int replicaExchangeUpdate(int commSize, int myRank, MPI_Comm comm, Aggregate &myAggregate, double temperature_array[])
{
	int exchangeAccepted = 0;
	
	// MASTER PROCESS 
	if(myRank == 0)
	{
		// Set-up exchange arrays to hold indices of send/receive processes
		int exchangeTemporary;
		int *exchangeReceieve	= new int[commSize];
    	int *exchangeSend		= new int[commSize];
		
		for (int i = 0; i < commSize; i ++)
		{
			exchangeReceieve[i] = i;
			exchangeSend[i]    	= i;
		}

		// Gather the current pre-exchange energies from all processes
		double localEnergies[commSize];

		for (int i = 1; i < commSize; i++) MPI_Recv(&localEnergies[i], 1, MPI_DOUBLE,i, 0, comm, MPI_STATUS_IGNORE);

		// Determine whether an exchange will occur between two processes 
		double expFactor,randFactor;

		for (int i = 1; i < (commSize-1); i ++)
	    {
	       	expFactor  = MIN(1.0, exp((localEnergies[exchangeReceieve[i+1]] - localEnergies[exchangeReceieve[i]])*((1/(temperature_array[i+1]))-(1/(temperature_array[i])))));
	       	randFactor = genrandBasic();

	       	if (randFactor <= expFactor)
	       	{
	       		exchangeTemporary    	= exchangeReceieve[i];
	       		exchangeReceieve[i]   	= exchangeReceieve[i+1];
	       		exchangeReceieve[i+1] 	= exchangeTemporary;
	       	}	
	    }

	    // Determine the send-to array 
		for (int i = 1; i < commSize; i++)
		{
			for(int j = 1; j < commSize; j++)
			{
				if(exchangeReceieve[j] == i) exchangeSend[i] = j;	
			}
		}
	       			
	    // Send both receive/send exchange arrays to the appropriate processes 
	    for (int i = 1; i < commSize; i++)
	    {
	       	MPI_Send(&exchangeSend[i], 1, MPI_INT, i, 0, comm);
	       	MPI_Send(&exchangeReceieve[i], 1, MPI_INT, i, 1, comm);
	    }

	    // Clean up
	    delete [] exchangeReceieve;
	    delete [] exchangeSend;	   
	}

	// SLAVE PROCESSES
	else
	{
		// Send local energies to Master 
		double localEnergy = myAggregate.getTotalEnergy();
		MPI_Send(&localEnergy, 1, MPI_DOUBLE, 0, 0, comm); 
		       
		// Find out if Replica Exchange Occurs
		int sendTo,receiveFrom;
		MPI_Recv(&sendTo, 1, MPI_INT, 0, 0, comm, MPI_STATUS_IGNORE);
		MPI_Recv(&receiveFrom, 1, MPI_INT, 0, 1, comm, MPI_STATUS_IGNORE);
		       		
		// Carry out the replica update by exchanging the coordinates of all the monomers in the aggregate
		// Note: Energies and other relevant quantities must be recalculated localy after the exchange   
		if (myRank != receiveFrom)
		{
			exchangeAccepted = 1;
			// Store aggregate coordinates in an array
			int coordinateArraySize 		= AGGREGATE_SIZE*POLYMER_LENGTH*3;
			double *coordinateArraySend 	= new double[coordinateArraySize];
			double *coordinateArrayReceive 	= new double[coordinateArraySize];

			for (int i = 0; i < AGGREGATE_SIZE; i++)
			{
				for(int j = 0; j < POLYMER_LENGTH; j++)
				{
					coordinateArraySend[(i*POLYMER_LENGTH*3)+(3*j)]		= (myAggregate.PolymerArray[i]-> MonomerArray[j]).getX();
					coordinateArraySend[(i*POLYMER_LENGTH*3)+(3*j)+1]	= (myAggregate.PolymerArray[i]-> MonomerArray[j]).getY();
					coordinateArraySend[(i*POLYMER_LENGTH*3)+(3*j)+2]	= (myAggregate.PolymerArray[i]-> MonomerArray[j]).getZ();
				}
			}

			// Exchange the coordinate arrays
	      	MPI_Sendrecv(coordinateArraySend, coordinateArraySize, MPI_DOUBLE, sendTo, 0, coordinateArrayReceive, coordinateArraySize, MPI_DOUBLE, receiveFrom, 0, comm, MPI_STATUS_IGNORE);

	      	// Update the positions of all monomers
	      	for (int i = 0; i < AGGREGATE_SIZE; i++)
			{
				for(int j = 0; j < POLYMER_LENGTH; j++)
				{
					(myAggregate.PolymerArray[i]-> MonomerArray[j]).setPosition(coordinateArrayReceive[(i*POLYMER_LENGTH*3)+(3*j)],coordinateArrayReceive[(i*POLYMER_LENGTH*3)+(3*j)+1],coordinateArrayReceive[(i*POLYMER_LENGTH*3)+(3*j)+2]);
				}
			}

			// Update the energies of the individual polymers as well as the interaction energy of the cluster
			myAggregate.recalculateAggregateEnergy();

			// Clean up
			delete [] coordinateArraySend;
			delete [] coordinateArrayReceive;
		}
	}
	return exchangeAccepted;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                               // OPTIMIZE DISPLACEMENT VALUES //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double optimizeDisplacementValues(int myRank, Aggregate &myAggregate, double myTemperature)
{
	int trainingComplete				= 0;
	double metropolisAcceptanceRatio 	= 0.0;
	double trainedDisplacement 			= SINGLE_DISPLACEMENT;

	while(trainingComplete == 0)
	{
		for (int i = 0; i < TRAINING_PHASE; i++) metropolisAcceptanceRatio += myAggregate.performSingleMetropolisUpdate(SPHERICAL_CONSTRAINT,trainedDisplacement,myTemperature,MAX_ENERGY);
		
		metropolisAcceptanceRatio = metropolisAcceptanceRatio/TRAINING_PHASE;
		
		if ((metropolisAcceptanceRatio < 0.6 && metropolisAcceptanceRatio > 0.5) || trainedDisplacement > 1.0 || trainedDisplacement < 0.01) trainingComplete = 1;
		else
		{
			if (metropolisAcceptanceRatio <= 0.5) trainedDisplacement = trainedDisplacement*0.9;
			else trainedDisplacement = trainedDisplacement/0.9;
			metropolisAcceptanceRatio  = 0.0;
		}
	} 
	
	return trainedDisplacement;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	                               // HISTOGRAM PRINTING //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void printEnergyHistograms(int myRank, ExtendedHistogram &myHistogram, double myTemperature, int fileIndex)
{

	// Determine histogram parameters
	double minEnergy,binWidth;
	int numberOfBins;

	numberOfBins 	= myHistogram.getNumberOfBins();
	minEnergy 		= myHistogram.getMin();
	binWidth 		= myHistogram.getBinWidth();

	// Open a write-only file
	MPI_File histoFile;
	char   mpiFileName[128];
    sprintf(mpiFileName, "canonicalHistograms_%d",fileIndex);
	MPI_File_open(MPI_COMM_WORLD, mpiFileName, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &histoFile);
	MPI_File_set_view(histoFile, myRank * (numberOfBins+1) * sizeof(double), MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL); 

	// Rank zero process prints out the energy ranges of each histogram bin
	if(myRank == 0)
	{
		double energyBin;
		double placeHolder = 0.0;
		MPI_File_write(histoFile ,&placeHolder, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);

		for(int i = 0; i < numberOfBins; i ++)
		{
			energyBin = minEnergy + (numberOfBins * i);
			MPI_File_write(histoFile ,&energyBin, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
		} 
	}
	else
	{
		double histogramValue;
		myHistogram.normalizeHistogram();
		MPI_File_write(histoFile ,&myTemperature, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);

		for(int i = 0; i < numberOfBins; i ++)
		{
			histogramValue = myHistogram.getValueAtIndex(i);
			MPI_File_write(histoFile ,&histogramValue, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
		}
	}
	MPI_File_close(&histoFile);
}

