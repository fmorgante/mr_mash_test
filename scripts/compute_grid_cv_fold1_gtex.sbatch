#!/bin/bash

#SBATCH --job-name=compute_grid_cv_fold1_gtex
#SBATCH --output=compute_grid_cv_fold1_gtex.%j.out
#SBATCH --error=compute_grid_cv_fold1_gtex.%j.err
#SBATCH --time=27:00:00
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --account=pi-mstephens

###Activate environment and set variables
source /scratch/midway3/fmorgante/miniconda3/etc/profile.d/conda.sh
conda activate py38dsc

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

###Cd into working directory
cd ../output/gtex_mr_mash_analysis/grid

###Extract effects
sos run ../../../../bioworkflows/multivariate-prediction/prior_grid.ipynb compute_grid \
                --analysis-units ../data_driven_matrices/analysis_units.txt \
                --datadir ../summary_stats --name fold_1 --table_name fold_1 \
                --bhat Bhat --sbhat Shat --expected-ncondition 49 \
                -c ../../../scripts/midway3_sos.yml -q midway3

###Deactivate environment
conda deactivate
