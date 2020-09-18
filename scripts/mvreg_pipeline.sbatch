#!/bin/bash

#SBATCH --job-name=mvreg_pipeline
#SBATCH --output=mvreg_pipeline.%j.out
#SBATCH --error=mvreg_pipeline.%j.err
#SBATCH --time=02:00:00
#SBATCH --partition=mstephens
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G

###Activate the environment
source activate py37dsc

###Set variables
dsc_outname=dsc_test
dsc_outdir="../output/${dsc_outname}"
inter_outdir="../output/${dsc_outname}_inter"
n_ind=90 ###needs to match the value in dsc
total_genes=100
pred_genes=10

###Load necessary version of R
module load gcc/6.2 texlive/2017 java/1.8 curl/7.50 mkl/2019.up1
export MKL_NUM_THREADS=1
export RHOME=/project2/mstephens/software/R-3.5.3-mkl
export PATH=$RHOME/bin:$PATH
export LD_LIBRARY_PATH=$RHOME/lib64/R/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH=$RHOME/lib64/R/lib:$LIBRARY_PATH

###Make some necessary directories
mkdir -p "${inter_outdir}/prior/sumstats"
mkdir -p "${inter_outdir}/prior/matrices"
mkdir -p "${inter_outdir}/misc"

###Run dsc to simulate genes and process data for data-driven covariance estimation
cd ../dsc
Rscript ../code/sample_testset_indeces.R --output "${inter_outdir}/misc/testset_indeces.rds" \
                                         --n "${n_ind}" \
                                         --prop_testset 0.2
                                         
./mvreg_all_genes_prior_09_11_20.dsc -o "${dsc_outdir}" --target sim_proc --replicate "${total_genes}" --host midway2.yml

###Prepare summary stats for data-driven covariance estimation
Rscript ../code/prepare_sumstats_for_ED_prior.R --input_path "${dsc_outdir}/univ_sumstats" \
                                                --output "${inter_outdir}/prior/sumstats/${dsc_outname}_effects_for_ED_prior.rds" \
                                                --n_random 4
                                                
###Run GTEx pipeline to compute the prior covariance matrices
sos run /project2/mstephens/fmorgante/gtexresults/workflows/mashr_flashr_workflow.ipynb prior \
        -c ../scripts/midway2_sos.yml \
        -q midway2 \
        --cwd "${inter_outdir}/prior/matrices" \
        --data "${inter_outdir}/prior/sumstats/${dsc_outname}_effects_for_ED_prior.rds" \
        --output-prefix "${dsc_outname}"

###Run dsc to simulate genes and process data for data-driven covariance estimation
./mvreg_all_genes_prior_09_11_20.dsc -o "${dsc_outdir}" --target fit_pred_score -s existing --replicate "${pred_genes}" --host midway2.yml

cd ../scripts

###Deactivate the environment
source deactivate



