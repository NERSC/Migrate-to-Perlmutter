#!/bin/bash
#SBATCH -q debug
#SBATCH -N 2
#SBATCH -C haswell
#SBATCH -t 00:10:00
#SBATCH -J my_job

srun -n 64 -c 2 --cpu-bind=cores ./mpi-hello
