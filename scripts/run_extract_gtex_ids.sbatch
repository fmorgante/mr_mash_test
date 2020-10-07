#!/bin/bash

#SBATCH --job-name=extract_gtex_ids
#SBATCH --output=extract_gtex_ids.%j.out
#SBATCH --error=extract_gtex_ids.%j.err
#SBATCH --time=00:01:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --account=pi-mstephens

###Load modules
module load R/3.5.1

###Run script
Rscript ../code/extract_gtex_ids.R
