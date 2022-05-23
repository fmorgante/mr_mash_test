#!/bin/bash

#SBATCH --job-name=prepare_sims_intermediate_files
#SBATCH --output=prepare_sims_intermediate_files.%j.out
#SBATCH --error=prepare_sims_intermediate_files.%j.err
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
declare -a arr=("mvreg_all_genes_prior_GTExrealX_indepV_B_1causalrespr10_indepReps" "mvreg_all_genes_prior_GTExrealX_indepV_indepB_r10_indepReps" "mvreg_all_genes_prior_GTExrealX_indepV_sharedB_2blocksr10_diffPVE_indepReps" "mvreg_all_genes_prior_GTExrealX_indepV_sharedB_3causalrespr10_indepReps" "mvreg_all_genes_prior_GTExrealX_indepV_sharedB_r10_indepReps")

for elem in "${arr[@]}"
do
  ###Decompress results files
  tar -xkf ../../output/"$elem"_inter.tar --directory ../../output/
  tar -xkf ../../output/"$elem".tar --directory ../../output/

  ###Process the results
  Rscript ../../code/paper_figures/sims_mvreg_all_genes_prior_GTExrealX_indepV_prep.R --input "$elem"

  ###Delete decompressed directories
  rm -rf ../../output/"$elem"_inter
  rm -rf ../../output/"$elem"
done

###Deactivate environment
conda deactivate
