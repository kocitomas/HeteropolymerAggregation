OBJS	= Histogram.o HistogramExtended.o MersenneTwister.o RandomNumberGenerator.o Vector.o Monomer.o Polymer.o Aggregate.o HelperFunctions.o
CC		= g++
MP		= mpic++
CFLAGS	= -c

## Link

all: test

AggregationHeteropolymer: $(OBJS) AggregationHeteropolymer.o
	$(MP) $(OBJS) AggregationHeteropolymer.o -o AggregationHeteropolymer



## Compile
MersenneTwister.o: MersenneTwister.h MersenneTwister.cpp
	$(CC) $(CFLAGS) MersenneTwister.cpp

RandomNumberGenerator.o: MersenneTwister.h RandomNumberGenerator.h RandomNumberGenerator.cpp
	$(CC) $(CFLAGS) RandomNumberGenerator.cpp

Vector.o: Vector.h Vector.cpp
	$(CC) $(CFLAGS) Vector.cpp

Monomer.o: Vector.h SimulationParameters.h Monomer.h Monomer.cpp
	$(CC) $(CFLAGS) Monomer.cpp
	
Polymer.o: Monomer.h RandomNumberGenerator.h Polymer.h Polymer.cpp
	$(CC) $(CFLAGS) Polymer.cpp

Aggregate.o: Polymer.h HistogramExtended.o Aggregate.h Aggregate.cpp
	$(CC) $(CFLAGS) Aggregate.cpp

Histogram.o: Histogram.h Histogram.cpp
	$(CC) $(CFLAGS) Histogram.cpp

HistogramExtended.o: Histogram.h HistogramExtended.h HistogramExtended.cpp
	$(CC) $(CFLAGS) HistogramExtended.cpp

HelperFunctions.o: HistogramExtended.h Aggregate.h HelperFunctions.h HelperFunctions.cpp
	$(MP) $(CFLAGS) HelperFunctions.cpp	

AggregationHeteropolymer.o: Aggregate.h HelperFunctions.h AggregationHeteropolymer.cpp
	$(MP) $(CFLAGS) AggregationHeteropolymer.cpp


clean:
	\rm *.o
