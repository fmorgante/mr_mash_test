#!/bin/bash

#SBATCH --job-name=mvreg_mrmash_vs_mtlasso_vs_utmost
#SBATCH --output=mvreg_mrmash_vs_mtlasso_vs_utmost.%j.out
#SBATCH --error=mvreg_mrmash_vs_mtlasso_vs_utmost.%j.err
#SBATCH --time=08:00:00
#SBATCH --partition=caslake
#SBATCH --account=pi-mstephens
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G

###Activate the environment
source /scratch/midway3/fmorgante/miniconda3/etc/profile.d/conda.sh
conda activate py38dsc

###Set variables
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
reps=20

###Run dsc
mkdir -p ../output/mvreg_mrmash_vs_mtlasso_vs_utmost
cd ../dsc

./mvreg_mrmash_vs_mtlasso_vs_utmost.dsc -o ../output/mvreg_mrmash_vs_mtlasso_vs_utmost --replicate ${reps} -c 2 --host midway3.yml

###Deactivate the environment
conda deactivate
