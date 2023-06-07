#!/bin/bash

#SBATCH --job-name=fit_mrmash_wtccc_sim
#SBATCH --output=fit_mrmash_wtccc_sim.%j.out
#SBATCH --error=fit_mrmash_wtccc_sim.%j.err
#SBATCH --time=36:00:00
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --account=pi-mstephens

###Activate environment and set variables
source /scratch/midway3/fmorgante/miniconda3/etc/profile.d/conda.sh
conda activate py38dsc

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

###Run R script
Rscript ../code/fit_mrmash_wtccc_sim.R

###Deactivate environment
conda deactivate
