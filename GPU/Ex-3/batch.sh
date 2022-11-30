#!/bin/bash
#SBATCH -N 1
#SBATCH -C gpu
#SBATCH -G 4
#SBATCH -t 5
# replace xxx below with the NERSC repo you plan to use, eg m1234
#SBATCH -A xxx
# uncomment the following lines to use the training reservation (Dec 1 only):
##SBATCH --reservation=pm_gpu_dec1
##SBATCH -A ntrain2

srun -n4 ./vec_add

