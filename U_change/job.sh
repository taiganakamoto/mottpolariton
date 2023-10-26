#!/bin/bash
#PBS -N sample-job
#PBS -q normal
#PBS -l
#PBS -j oe

cd $PBS_O_WORKDIR

export JULIA_NUM_THREADS=28

julia ./opcon.jl
