#!/bin/bash
#PBS -N sample-job
#PBS -q normal
#PBS -l nodes=qcrt04.ad12.riken.jp:ppn=28 
#PBS -j oe

cd $PBS_O_WORKDIR

export JULIA_NUM_THREADS=28

julia ./gene.jl 4
