#!/bin/bash
#SBATCH -q debug
#SBATCH -N 2
#SBATCH -C cpu
#SBATCH -t 00:10:00
#SBATCH -J run-xthi

export OMP_NUM_THREADS=8
export OMP_PROC_BIND=spread
export OMP_PLACES=threads

echo "== run 1 ===="
srun -n 8 -c 64 --cpu-bind=cores ./xthi

echo "== run 2, sorted ===="
srun -n 8 -c 64 --cpu-bind=cores ./xthi |sort -k4n,6n


echo "====run 3, with display_affinity============="
export OMP_DISPLAY_AFFINITY=true
export OMP_AFFINITY_FORMAT="host=%H, pid=%P, thread_num=%n, thread affinity=%A"
srun -n 8 -c 64 --cpu-bind=cores ./xthi |sort -k4n,6n
