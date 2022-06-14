#!/bin/bash

#SBATCH --job-name=mrmash_enet_mtlasso_score_10genes_gtex
#SBATCH --output=mrmash_enet_mtlasso_score_10genes_gtex.%j.out
#SBATCH --error=mrmash_enet_mtlasso_score_10genes_gtex.%j.err
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

###Run SoS
sos run /project/mstephens/fmorgante/mr_mash_test/code/prediction_score.ipynb \
    --analysis-units /project/mstephens/fmorgante/mr_mash_test/data/gtex-v8-manifest-2ormore-tissues-nopath-nosuffix-mtlasso.txt \
    --data-dir /project/mstephens/fmorgante/mr_mash_test/data/cis_eqtl_analysis_ready \
    --mrmash-model /project/mstephens/fmorgante/mr_mash_test/output/gtex_mr_mash_analysis/prediction \
    --mtlasso-model /project/mstephens/fmorgante/mr_mash_test/output/gtex_mr_mash_analysis/mtlasso \
    --analysis-stage first_pass \
    --data-suffix GTEx_V8 \
    --wd /project/mstephens/fmorgante/mr_mash_test/output/gtex_mr_mash_analysis/prediction_score_mtlasso \
    -c midway3_sos.yml -q midway3

###Deactivate environment
conda deactivate

