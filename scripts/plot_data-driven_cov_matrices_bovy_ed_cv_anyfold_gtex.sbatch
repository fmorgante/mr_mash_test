#!/bin/bash

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

###Set variables
fold_num="${1}"

###Cd into working directory
mkdir -p ../output/gtex_mr_mash_analysis/data_driven_matrices/plots
cd ../output/gtex_mr_mash_analysis/data_driven_matrices/plots

###Denoise estimated matrices by UD (ED method)
sos run ../../../../../bioworkflows/multivariate-fine-mapping/mixture_prior.ipynb plot_U \
		--model-data "../output/fold_${fold_num}.ed_bovy.rds" --name "" --cwd . \
		--remove-label False

###Deactivate environment
conda deactivate
