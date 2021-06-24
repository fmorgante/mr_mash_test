#!/bin/bash

#SBATCH --job-name=compute_data-driven_cov_matrices_ud_ed_cv_fold2_gtex
#SBATCH --output=compute_data-driven_cov_matrices_ud_ed_cv_fold2_gtex.%j.out
#SBATCH --error=compute_data-driven_cov_matrices_ud_ed_cv_fold2_gtex.%j.err
#SBATCH --time=36:00:00
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

###Denoise estimated matrices by UD (ED method)
sos run ../../../../bioworkflows/multivariate-fine-mapping/mixture_prior.ipynb ud \
                --datadir ../summary_stats --name fold_2 --ud-tol-lik 1e-3 \
		--mixture-components flash flash_nonneg pca \
                -c ../../../scripts/midway3_sos.yml -q midway3

###Deactivate environment
conda deactivate
