#!/bin/bash

#SBATCH --job-name=compute_sumstats_cv_gtex
#SBATCH --output=compute_sumstats_cv_gtex.%j.out
#SBATCH --error=compute_sumstats_cv_gtex.%j.err
#SBATCH --time=12:00:00
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
sos run /project/mstephens/fmorgante/bioworkflows/misc/cv_preprocessing.ipynb \
    --analysis-units ../data/gtex-v8-manifest-nopath-nosuffix-10genes.txt \
    --sample-partition ../data/gtex-v8-ids-folds.txt \
    --data-dir ../data/cis_eqtl_analysis_ready \
    --data-suffix GTEx_V8.rds \
    --name sumstats_cv \
    --wd ../output/gtex_mr_mash_analysis/summary_stats \
    --nthreads 1 \
    -c midway3_sos.yml -q midway3
    
###Deactivate environment
conda deactivate

