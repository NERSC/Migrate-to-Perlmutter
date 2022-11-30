#!/bin/bash
#SBATCH -N 1
#SBATCH -C gpu
#SBATCH -G 4
#SBATCH -t 5
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=2
# replace xxx below with the NERSC repo you plan to use, eg m1234
#SBATCH -A xxx
# uncomment the following lines to use the training reservation (Dec 1 only):
##SBATCH --reservation=pm_gpu_dec1
##SBATCH -A ntrain2

lscpu | grep NUMA
echo -e "\n"
lstopo | grep -e "NUMANode" -e "3D"
echo -e "\n"
srun -n8 --cpu-bind=cores --gpu-bind=closest ./vec_add

