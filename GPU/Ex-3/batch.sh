#!/bin/bash
#SBATCH -N 1
#SBATCH -C gpu
#SBATCH -G 4
#SBATCH -t 5
#SBATCH -A xxx
#SBATCH --reservation=xxxx

srun -n4 ./vec_add