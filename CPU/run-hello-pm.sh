#!/bin/bash
#SBATCH -q debug
#SBATCH -N 2
#SBATCH -C cpu
#SBATCH -t 00:10:00
#SBATCH -J my_job

srun -n 64 -c 8 --cpu-bind=cores ./mpi-hello
