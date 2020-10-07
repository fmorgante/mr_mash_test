#!/bin/bash

#SBATCH --job-name=get_full_X_gtex_v8_manifest
#SBATCH --output=get_full_X_gtex_v8_manifest.%j.out
#SBATCH --error=get_full_X_gtex_v8_manifest.%j.err
#SBATCH --time=05:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --account=pi-mstephens

###Load modules
module load R/3.5.1

###Run script
Rscript ../code/get_full_X_gtex_v8_manifest.R --n_cores 4
