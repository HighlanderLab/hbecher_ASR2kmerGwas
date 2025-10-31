#!/bin/bash
#SBATCH --job-name=uni-ST
#SBATCH --output=slurm-%A.out
#SBATCH --error=slurm-%A.err
#SBATCH -A k10170
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=workq
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH --time=12:00:00
#SBATCH --exclusive

##OpenMP settings:
#export OMP_NUM_THREADS=1


#cat run_gwas.sh

#run the application:

#source /home/omairyrm/miniconda3/etc/profile.d/conda.sh
#conda activate r-env
source /scratch/omairyrm/miniconda3/etc/profile.d/conda.sh
conda init
conda activate
conda activate r-env

Rscript 00RUNME.R



