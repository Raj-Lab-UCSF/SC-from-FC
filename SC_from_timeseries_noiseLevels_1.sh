#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#### Specify job name
#$ -N SC_from_timeseries
#### Output file
#$ -o /wynton/home/rajlab/fabdelnour/jobs/SC_from_fitSGM/JobOutput/$JOB_NAME_$JOB_ID.out
#### Error file
#$ -e /wynton/home/rajlab/fabdelnour/jobs/SC_from_fitSGM/JobOutput/$JOB_NAME_$JOB_ID.out
#### number of cores
#$ -pe smp 50
#### memory per core
#$ -l mem_free=2G
#### Maximum run time
#$ -l h_rt=150:00:00

module load matlab
matlab -batch "runjob_noiseLevels_1"