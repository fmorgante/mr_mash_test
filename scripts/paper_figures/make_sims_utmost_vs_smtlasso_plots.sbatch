#!/bin/bash

#SBATCH --job-name=make_sims_utmost_vs_smtlasso_plots
#SBATCH --output=make_sims_utmost_vs_smtlasso_plots.%j.out
#SBATCH --error=make_sims_utmost_vs_smtlasso_plots.%j.err
#SBATCH --time=00:30:00
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
  
###Process the results
Rscript ../../code/paper_figures/sims_mvreg_utmost_vs_smtlasso_plot.R

###Deactivate environment
conda deactivate
