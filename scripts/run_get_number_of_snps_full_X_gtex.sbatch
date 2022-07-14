#!/bin/bash

#SBATCH --job-name=get_number_of_snps_full_X_gtex
#SBATCH --output=get_number_of_snps_full_X_gtex.%j.out
#SBATCH --error=get_number_of_snps_full_X_gtex.%j.err
#SBATCH --time=05:00:00
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --account=pi-mstephens

###Activate environment and set variables
source /scratch/midway3/fmorgante/miniconda3/etc/profile.d/conda.sh
conda activate py38dsc

###Run script
Rscript ../code/get_number_of_snps_full_X_gtex.R --n_cores 4

###Deactivate environment
conda deactivate
