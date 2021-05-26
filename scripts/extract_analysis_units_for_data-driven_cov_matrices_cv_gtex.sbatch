#!/bin/bash

#SBATCH --job-name=extract_analysis_units_for_data-driven_cov_matrices_cv_gtex
#SBATCH --output=extract_analysis_units_for_data-driven_cov_matrices_cv_gtex.%j.out
#SBATCH --error=extract_analysis_units_for_data-driven_cov_matrices_cv_gtex.%j.err
#SBATCH --time=00:30:00
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --account=pi-mstephens

###Cd into working directory
cd ../output/gtex_mr_mash_analysis/data_driven_matrices

###Prepare a list of genes to be used in the correct format
find ../summary_stats/ -type f -name '*.rds' -exec basename {} .rds \; > analysis_units.txt
