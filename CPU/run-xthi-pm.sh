#!/bin/bash
#SBATCH -q debug
#SBATCH -N 2
#SBATCH -C cpu
#SBATCH -t 00:10:00
#SBATCH -J run-xthi

export OMP_NUM_THREADS=8
export OMP_PROC_BIND=spread
export OMP_PLACES=threads

srun -n 8 -c 64 --cpu-bind=cores ./xthi
