#!/bin/bash -l

#$ -S /bin/bash

### you can name your job, otherwise it will be listed by the name of this script ###
#$ -N AggregateStudyHeteropolymer

### Tell queue system to use the directory of the submit script ###
### as the current working directory ###
#$ -cwd

### set the proper queue for your job, either smsgroup1 or smsgroup2 ###
#$ -q smsgroup3

### set the parallel environment for the queue (g1mpi or g2mpi) and request processors ###
#$ -pe g3mpi 92
### for the best performance, mpi jobs should run from the infiniband scratch dir ###
### your scratch dir is set up under the environment variable $IBHOME ###
### create directory to run the job via infiniband ###

#RUNDIR=${IBHOME}/AggregateStudy-$JOB_ID
#mkdir $RUNDIR 

### create a directory for the results in your home directory to make them available ###
### to your dekstop machine ###

RESULTDIR=`pwd`/AggregateStudy-$JOB_ID
RUNDIR=${RESULTDIR}
mkdir $RESULTDIR

### copy code and any needed input files to RUNDIR ###
cp temperatures SimulationParameters.h AggregationHeteropolymer $RUNDIR

### change current working directory to job directory ###
cd $RUNDIR 

### execute binary ###
time mpirun -np $NSLOTS --mca btl sm --mca btl tcp,self --mca btl_tcp_if_include em1  AggregationHeteropolymer

###copy files back to $RESULTDIR in home directory ###
#cp * $RESULTDIR 

### cleanup $IBHOME since space is limited ###
### I've commented this out for now to make sure everything works properly before stuff is deleted. Can be manually deleted later ###
#rm -rf $RUNDIR
