#!/bin/bash

#SBATCH --job-name=mrmash_joint_weights_update_fold1_gtex
#SBATCH --output=mrmash_joint_weights_update_fold1_gtex.%j.out
#SBATCH --error=mrmash_joint_weights_update_fold1_gtex.%j.err
#SBATCH --time=07:00:00
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

cd ../output/gtex_mr_mash_analysis/weights

###Create analysis units file
find ../prediction/fold_1 -maxdepth 1 -type f -name '*.rds' -exec basename {} \; > tmp_fold_1.txt
cut -d"." -f1-2 tmp_fold_1.txt > analysis_units_weights_fold_1.txt
rm tmp_fold_1.txt

###Run SoS
sos run /project/mstephens/fmorgante/bioworkflows/multivariate-prediction/mrmash.ipynb joint_weights_update \
    --analysis-units analysis_units_weights_fold_1.txt \
    --data-dir ../prediction/fold_1  \
    --data-suffix GTEx_V8_fold_1_mrmash_first_pass.rds \
    --name fold_1 \
    --wd . \
    -c ../../../scripts/midway3_sos.yml -q midway3

###Deactivate environment
conda deactivate
