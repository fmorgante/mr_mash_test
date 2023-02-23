#!/bin/bash

#SBATCH --job-name=prepare_sims_utmost_vs_smtlasso_intermediate_files
#SBATCH --output=prepare_sims_utmost_vs_smtlasso_intermediate_files.%j.out
#SBATCH --error=prepare_sims_utmost_vs_smtlasso_intermediate_files.%j.err
#SBATCH --time=04:00:00
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

###Loop over scenarios
declare -a arr=("indepX_indepV_sharedB_allr_norm" "indepX_indepV_indepB_allr_norm" "indepX_indepV_B_1causalr" "indepX_indepV_sharedB_3causalrespr10" "indepX_indepV_sharedB_2blocksr_diffPVE")

for elem in "${arr[@]}"
do
  ###Process the results
  Rscript ../../code/paper_figures/sims_mvreg_utmost_vs_smtlasso_prep.R --input "$elem"
done

###Deactivate environment
conda deactivate
