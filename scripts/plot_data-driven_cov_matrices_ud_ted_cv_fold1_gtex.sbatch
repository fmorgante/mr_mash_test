#!/bin/bash

#SBATCH --job-name=plot_data-driven_cov_matrices_ud_ted_cv_fold1_gtex
#SBATCH --output=plot_data-driven_cov_matrices_ud_ted_cv_fold1_gtex.%j.out
#SBATCH --error=plot_data-driven_cov_matrices_ud_ted_cv_fold1_gtex.%j.err
#SBATCH --time=01:00:00
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
mkdir -p ../output/gtex_mr_mash_analysis/data_driven_matrices/plots
cd ../output/gtex_mr_mash_analysis/data_driven_matrices/plots

###Denoise estimated matrices by UD (ED method)
sos run ../../../../../bioworkflows/multivariate-fine-mapping/mixture_prior.ipynb plot_U \
		--model-data ../output/fold_1.ted_unconstrained.rds --name "" --cwd . \
		--remove-label False

###Deactivate environment
conda deactivate
