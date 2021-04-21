#!/bin/bash

#SBATCH --job-name=get_missing_Y_rate_by_tissue_gtex_v8_manifest
#SBATCH --output=get_missing_Y_rate_by_tissue_gtex_v8_manifest.%j.out
#SBATCH --error=get_missing_Y_rate_by_tissue_gtex_v8_manifest.%j.err
#SBATCH --time=05:00:00
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G
#SBATCH --account=pi-mstephens

###Activate conda environment
source /scratch/midway3/fmorgante/miniconda3/etc/profile.d/conda.sh
conda activate py38dsc

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

###Run script
Rscript ../code/get_missing_Y_rate_by_tissue_gtex_v8_manifest.R --n_cores 4

###Deactivate conda environment
conda deactivate
