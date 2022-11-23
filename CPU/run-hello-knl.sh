#!/bin/bash
#SBATCH -q debug
#SBATCH -N 1
#SBATCH -C haswell
#SBATCH -t 00:10:00
#SBATCH -J my_job

echo "from serial-hello"
./serial-hello

echo "from mpi-hello"
srun -n 8 -c 8 --cpu-bind=cores ./mpi-hello
