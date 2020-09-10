#!/bin/bash

#SBATCH --job-name=mvreg_pipeline
#SBATCH --output=mvreg_pipeline.%j.out
#SBATCH --error=mvreg_pipeline.%j.err
#SBATCH --time=01:00:00
#SBATCH --partition=mstephens
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G

###Activate the environment
source activate py37dsc

###Set variables
dsc_outname=dsc_test
dsc_outdir="../output/${dsc_outname}"
prior_outdir="../output/${dsc_outname}_prior"

###Load necessary version of R
module load gcc/6.2 texlive/2017 java/1.8 curl/7.50 mkl/2019.up1
export MKL_NUM_THREADS=1
export RHOME=/project2/mstephens/software/R-3.5.3-mkl
export PATH=$RHOME/bin:$PATH
export LD_LIBRARY_PATH=$RHOME/lib64/R/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH=$RHOME/lib64/R/lib:$LIBRARY_PATH

###Run dsc to simulate genes
cd ../dsc
./simulate_process_data_for_mvreg_09_09_20.dsc -o "${dsc_outdir}" -c 4 

###Prepare summary stats for data-driven covariance estimation
mkdir -p "${prior_outdir}/sumstats"
mkdir -p "${prior_outdir}/matrices"
Rscript ../code/prepare_sumstats_for_ED_prior.R --input_path "${dsc_outdir}/univ_sumstats" \
                                                --output "${prior_outdir}/sumstats/effects_for_ED_prior.rds" \
                                                --n_random 4
                                                
###Run GTEx pipeline to compute the prior covariance matrices
sos run /project2/mstephens/fmorgante/gtexresults/workflows/mashr_flashr_workflow.ipynb prior \
        -j 4 \
        --cwd "${prior_outdir}/matrices" \
        --data "${prior_outdir}/sumstats/effects_for_ED_prior.rds"
                                        
###Deactivate the environment
source deactivate



