#!/bin/bash
#SBATCH --job-name=LPA
#SBATCH --account=co_noneq
#SBATCH --partition=savio2
#SBATCH --time=150:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --output="/global/scratch/users/seokjinmoon/logs/job_%j.out"

module load gcc
module load openmpi
module load parallel/20220522

echo $SLURM_CPUS_ON_NODE
parallel -j 5 < rate.lst
