#!/bin/bash

#SBATCH --job-name=prepare_sims_missing_Y_intermediate_files
#SBATCH --output=prepare_sims_missing_Y_intermediate_files.%j.out
#SBATCH --error=prepare_sims_missing_Y_intermediate_files.%j.err
#SBATCH --time=24:00:00
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
declare -a arr=("mvreg_all_genes_prior_GTExrealX_indepV_B_1causalrespr10_indepReps_missing_Y" "mvreg_all_genes_prior_GTExrealX_indepV_indepB_r10_indepReps_missing_Y" "mvreg_all_genes_prior_GTExrealX_indepV_sharedB_2blocksr10_diffPVE_indepReps_missing_Y" "mvreg_all_genes_prior_GTExrealX_indepV_sharedB_3causalrespr10_indepReps_missing_Y" "mvreg_all_genes_prior_GTExrealX_indepV_sharedB_r10_indepReps_missing_Y")

for elem in "${arr[@]}"
do
  ###Decompress results files
  tar -xf ../../output/"$elem"_inter.tar --directory /scratch/midway3/fmorgante/
  tar -xf ../../output/"$elem".tar --directory /scratch/midway3/fmorgante/

  ###Process the results
  Rscript ../../code/paper_figures/sims_mvreg_all_genes_prior_GTExrealX_indepV_missing_Y_prep.R --input /scratch/midway3/fmorgante/"$elem" --output ../../output/sims_paper_figures_inter/"$elem"

  ###Delete decompressed directories
  rm -rf /scratch/midway3/fmorgante/"$elem"_inter
  rm -rf /scratch/midway3/fmorgante/"$elem"
done

###Deactivate environment
conda deactivate
