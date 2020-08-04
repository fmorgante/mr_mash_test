#!/bin/bash

#SBATCH --job-name=dsc_mvreg_07_29_20
#SBATCH --output=dsc_mvreg_07_29_20.%j.out
#SBATCH --error=dsc_mvreg_07_29_20.%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=8G
#SBATCH --account=pi-yangili1

###Activate conda environment
source activate py37dsc

###Load appropriate R version (needed because dsc check for R dependencies)
module load gcc/6.2 texlive/2017 java/1.8 curl/7.50 mkl/2019.up1
export MKL_NUM_THREADS=1
export RHOME=/project2/mstephens/software/R-3.5.3-mkl
export PATH=$RHOME/bin:$PATH
export LD_LIBRARY_PATH=$RHOME/lib64/R/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH=$RHOME/lib64/R/lib:$LIBRARY_PATH

###Run dsc
./mvreg_07_29_20.dsc -o ../output/dsc_07_29_20 --host midway2.yml -c 5

###Deactivate conda environment
source deactivate
