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
sos run /project/mstephens/fmorgante/bioworkflows/multivariate-prediction/mrmash.ipynb mr_mash \
    --analysis-units "../data/gtex-v8-manifest-2ormore-tissues-nopath-nosuffix-batch${btc}.txt" \
    --data-dir ../data/cis_eqtl_analysis_ready \
    --data-suffix GTEx_V8.rds \
    --name "fold_3_batch_${btc}" \
    --wd ../output/gtex_mr_mash_analysis \
    --prior-grid ../output/gtex_mr_mash_analysis/grid/fold_3_grid.rds \
    --prior-matrices ../output/gtex_mr_mash_analysis/data_driven_matrices/output/fold_3.ted_unconstrained.rds \
    --sample-partition ../data/gtex-v8-ids-folds.txt \
    --fold 3 \
    --nthreads 2 \
    --update-w0 TRUE \
    --w0-threshold 1e-8 \
    --save-model TRUE \
    --glmnet-pred TRUE \
    -c midway3_sos.yml -q midway3

###Deactivate environment
conda deactivate
