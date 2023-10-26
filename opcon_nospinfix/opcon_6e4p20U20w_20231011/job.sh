#!/bin/bash
#PBS -N sample-job
#PBS -q normal
#PBS -l nodes=qcrt02.ad12.riken.jp:ppn=1
#PBS -j oe

cd $PBS_O_WORKDIR

export JULIA_NUM_THREADS=1

julia ./opcon.jl

