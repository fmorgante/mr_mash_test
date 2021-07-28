#!/bin/bash

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

btc=${1}

###Run SoS
sos run /project/mstephens/fmorgante/bioworkflows/multivariate-prediction/mrmash.ipynb joint_weights_update \
    --data-dir ../output/gtex_mr_mash_analysis/prediction/fold_4 \
    --data-suffix "GTEx_V8_fold_4_batch_${btc}_mrmash.first_pass.rds" \
    --name "fold_4_batch_${btc}" \
    --wd ../output/gtex_mr_mash_analysis \
    -c midway3_sos.yml -q midway3

###Deactivate environment
conda deactivate
