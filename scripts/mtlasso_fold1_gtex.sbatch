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

###Run SoS
sos run /project/mstephens/fmorgante/bioworkflows/multivariate-prediction/mtlasso.ipynb mtlasso \
    --analysis-units "../data/gtex-v8-manifest-2ormore-tissues-nopath-nosuffix-test.txt" \
    --data-dir ../data/cis_eqtl_analysis_ready \
    --data-suffix GTEx_V8.rds \
    --name "fold_1" \
    --wd ../output/gtex_mr_mash_analysis \
    --sample-partition ../data/gtex-v8-ids-folds.txt \
    --fold 1 \
    --mtlasso-script ../code/fit_mtlasso_missing_Y.py \
    --conda-env py38dsc \
    --grid-length 50 \
    -c midway3_sos.yml -q midway3

###Deactivate environment
conda deactivate
