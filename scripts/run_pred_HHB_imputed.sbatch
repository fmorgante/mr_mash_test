#!/bin/bash

#SBATCH --job-name=pred_HHB_imputed
#SBATCH --output=pred_HHB_imputed.%j.out
#SBATCH --error=pred_HHB_imputed.%j.err
#SBATCH --time=12:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --account=pi-mstephens

###Load modules
module load gcc/6.2 texlive/2017 java/1.8 curl/7.50 mkl/2019.up1
export MKL_NUM_THREADS=1
export RHOME=/project2/mstephens/software/R-3.5.3-mkl
export PATH=$RHOME/bin:$PATH
export LD_LIBRARY_PATH=$RHOME/lib64/R/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH=$RHOME/lib64/R/lib:$LIBRARY_PATH

###Run script
Rscript ../code/pred_HHB_imputed.R
