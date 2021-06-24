#!/bin/bash

#SBATCH --job-name=extract_effects_for_data-driven_cov_matrices_cv_fold4_gtex
#SBATCH --output=extract_effects_for_data-driven_cov_matrices_cv_fold4_gtex.%j.out
#SBATCH --error=extract_effects_for_data-driven_cov_matrices_cv_fold4_gtex.%j.err
#SBATCH --time=01:00:00
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --account=pi-mstephens

###Activate environment and set variables
source /scratch/midway3/fmorgante/miniconda3/etc/profile.d/conda.sh
conda activate py38dsc

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

###Cd into working directory
cd ../output/gtex_mr_mash_analysis/data_driven_matrices

###Extract effects
sos run ../../../../bioworkflows/multivariate-fine-mapping/mixture_prior.ipynb extract_effects \
                --analysis-units analysis_units.txt \
                --datadir ../summary_stats --name fold_4 --table_name fold_4 \
                --bhat Bhat --sbhat Shat --expected-ncondition 49 \
                -c ../../../scripts/midway3_sos.yml -q midway3

###Deactivate environment
conda deactivate
