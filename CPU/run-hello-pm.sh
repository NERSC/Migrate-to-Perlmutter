#!/bin/bash
#SBATCH -q debug
#SBATCH -N 2
#SBATCH -C cpu
#SBATCH -t 00:10:00
#SBATCH -J my_job
## uncomment the following lines to use the training reservation (Mar 10 only):
##SBATCH --reservation=pm_gpu_mar10
##SBATCH -A ntrain8

srun -n 64 -c 8 --cpu-bind=cores ./mpi-hello
