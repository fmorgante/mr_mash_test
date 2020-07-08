#!/bin/bash

#SBATCH --job-name=assess_mrmash_speed
#SBATCH --output=assess_mrmash_speed.%j.out
#SBATCH --error=assess_mrmash_speed.%j.err
#SBATCH --time=12:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --account=pi-yangili1

###Load necessary modules
module load gcc/6.2 texlive/2017 java/1.8 curl/7.50 mkl/2019.up1
export MKL_NUM_THREADS=1
export RHOME=/project2/mstephens/software/R-3.5.3-mkl
export PATH=$RHOME/bin:$PATH
export LD_LIBRARY_PATH=$RHOME/lib64/R/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH=$RHOME/lib64/R/lib:$LIBRARY_PATH

###Run R script
##No singletons
Rscript assess_mrmash_speed.R --nthreads 4 --standardize FALSE --w0_threshold 1e-8 --singletons FALSE
                     
echo -e "\n\n\nNO SINGLETON ANALYSIS FINISHED \n\n\n"

##With singletons
Rscript assess_mrmash_speed.R --nthreads 4 --standardize FALSE --w0_threshold 1e-8 --singletons TRUE

echo -e "\n\n\nSINGLETON ANALYSIS FINISHED\n\n\n"
