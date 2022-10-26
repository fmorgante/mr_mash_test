#!/bin/bash

#SBATCH --job-name=make_gtex_tPCA_heatmap_fold_1
#SBATCH --output=make_gtex_tPCA_heatmap_fold_1.%j.out
#SBATCH --error=make_gtex_tPCA_heatmap_fold_1.%j.err
#SBATCH --time=00:30:00
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
  
###Process the results
Rscript ../../code/paper_figures/gtex_tPCA_heatmap_fold_1.R

###Deactivate environment
conda deactivate
