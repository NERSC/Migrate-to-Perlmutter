#!/bin/bash
#SBATCH -N 1
#SBATCH -C gpu
#SBATCH -G 4
#SBATCH -t 5
#SBATCH -A xxx

module load PrgEnv-nvidia

make clean ; make

srun -n 2 ./bcast_from_device

