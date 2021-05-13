#!/bin/bash

#SBATCH --job-name=split_folds_gtex_ids
#SBATCH --output=split_folds_gtex_ids.%j.out
#SBATCH --error=split_folds_gtex_ids.%j.err
#SBATCH --time=00:01:00
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --account=pi-mstephens

###Activate environment and set variables
source /scratch/midway3/fmorgante/miniconda3/etc/profile.d/conda.sh
conda activate py38dsc

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

###Run script
Rscript ../code/split_folds_gtex_ids.R

###Deactivate environment 
conda deactivate
