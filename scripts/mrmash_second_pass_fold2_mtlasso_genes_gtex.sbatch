#!/bin/bash

#SBATCH --job-name=mrmash_second_pass_fold2_mtlasso_genes_gtex
#SBATCH --output=mrmash_second_pass_fold2_mtlasso_genes_gtex.%j.out
#SBATCH --error=mrmash_second_pass_fold2_mtlasso_genes_gtex.%j.err
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
sos run /project/mstephens/fmorgante/bioworkflows/multivariate-prediction/mrmash.ipynb mr_mash \
    --analysis-units ../data/gtex-v8-manifest-2ormore-tissues-nopath-nosuffix-mtlasso.txt \
    --data-dir ../data/cis_eqtl_analysis_ready  \
    --data-suffix GTEx_V8.rds \
    --name fold_2_batch_1 \
    --wd ../output/gtex_mr_mash_analysis \
    --prior-grid ../output/gtex_mr_mash_analysis/grid/fold_2_grid.rds \
    --prior-matrices ../output/gtex_mr_mash_analysis/data_driven_matrices/output/fold_2.ted_unconstrained.rds \
    --prior-weights ../output/gtex_mr_mash_analysis/weights/fold_2_batch_1_updated_weights.rds \
    --sample-partition ../data/gtex-v8-ids-folds.txt \
    --fold 2 \
    --nthreads 4 \
    --update-w0 FALSE \
    --w0-threshold 0 \
    --save-model TRUE \
    -c midway3_sos.yml -q midway3

###Deactivate environment
conda deactivate
