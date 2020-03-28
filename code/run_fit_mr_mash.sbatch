#!/bin/bash

#SBATCH --job-name=fit_mr_mash
#SBATCH --output=fit_mr_mash.%j.out
#SBATCH --error=fit_mr_mash.%j.err
#SBATCH --time=12:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --account=pi-yangili1

###Load necessary module
module load R/3.5.1

###Request cpus for openBLAS
export OPENBLAS_NUM_THREADS=4

###Run R script
##Independent predictors, shared effects
Rscript fit_mr_mash.R --n 600 --p 1000 --p_causal 50 --r 5 --pve 0.5 --sigma_offdiag 1 --sigma_scale 0.8 \
                      --gamma_offdiag 0 --gamma_scale 0.8 --V_offdiag 0.2 --V_scale 0 --prop_testset 0.2 \
                      --update_w0 TRUE --update_w0_method EM --standardize TRUE --verbose FALSE --update_V TRUE \
                      --outdir ./output/
                      
##Independent predictors, independent effects
Rscript fit_mr_mash.R --n 600 --p 1000 --p_causal 50 --r 5 --pve 0.5 --sigma_offdiag 0 --sigma_scale 0.8 \
                      --gamma_offdiag 0 --gamma_scale 0.8 --V_offdiag 0.2 --V_scale 0 --prop_testset 0.2 \
                      --update_w0 TRUE --update_w0_method EM --standardize TRUE --verbose FALSE --update_V TRUE \
                      --outdir ./output/

##Correlated predictors, shared effects
Rscript fit_mr_mash.R --n 600 --p 1000 --p_causal 50 --r 5 --pve 0.5 --sigma_offdiag 1 --sigma_scale 0.8 \
                      --gamma_offdiag 0.5 --gamma_scale 0.8 --V_offdiag 0.2 --V_scale 0 --prop_testset 0.2 \
                      --update_w0 TRUE --update_w0_method EM --standardize TRUE --verbose FALSE --update_V TRUE \
                      --outdir ./output/
                      
##Correlated predictors, independent effects
Rscript fit_mr_mash.R --n 600 --p 1000 --p_causal 50 --r 5 --pve 0.5 --sigma_offdiag 0 --sigma_scale 0.8 \
                      --gamma_offdiag 0.5 --gamma_scale 0.8 --V_offdiag 0.2 --V_scale 0 --prop_testset 0.2 \
                      --update_w0 TRUE --update_w0_method EM --standardize TRUE --verbose FALSE --update_V TRUE \
                      --outdir ./output/
