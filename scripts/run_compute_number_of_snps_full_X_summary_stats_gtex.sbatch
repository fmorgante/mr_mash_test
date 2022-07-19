#!/bin/bash

#SBATCH --job-name=compute_number_of_snps_full_X_summary_stats_gtex
#SBATCH --output=compute_number_of_snps_full_X_summary_stats_gtex.%j.out
#SBATCH --error=compute_number_of_snps_full_X_summary_stats_gtex.%j.err
#SBATCH --time=00:10:00
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --account=pi-mstephens

###Activate environment
source /scratch/midway3/fmorgante/miniconda3/etc/profile.d/conda.sh
conda activate py38dsc

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

###Run script
Rscript ../code/compute_number_of_snps_full_X_summary_stats_gtex.R

###Deactivate the environment
conda deactivate
