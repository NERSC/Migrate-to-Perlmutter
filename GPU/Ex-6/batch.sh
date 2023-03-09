#!/bin/bash
#SBATCH -N 1
#SBATCH -C gpu
#SBATCH -G 4
#SBATCH -t 5
## replace xxx below with the NERSC repo you plan to use, eg m1234
#SBATCH -A xxx
## uncomment the following lines to use the training reservation (Mar 10 only):
##SBATCH --reservation=pm_gpu_mar10
##SBATCH -A ntrain8

module load PrgEnv-nvidia

make clean ; make

srun -n 2 ./bcast_from_device

