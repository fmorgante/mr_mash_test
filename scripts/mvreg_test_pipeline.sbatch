#!/bin/bash

#SBATCH --job-name=mvreg_test_pipeline
#SBATCH --output=mvreg_test_pipeline.%j.out
#SBATCH --error=mvreg_test_pipeline.%j.err
#SBATCH --time=04:00:00
#SBATCH --partition=mstephens
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G

################
####Common steps
################

###Activate the environment
source activate py37dsc

###Set variables
dsc_script=./mvreg_test.dsc
dsc_outname=`echo "${dsc_script}" | cut -d'/' -f 2 | cut -d'.' -f 1`
dsc_outdir="../output/${dsc_outname}"
inter_outdir="../output/${dsc_outname}_inter"
n_testind=168
total_genes=10
pred_genes=1
n_reps=5


###Load necessary version of R
module load gcc/6.2 texlive/2017 java/1.8 curl/7.50 mkl/2019.up1
export MKL_NUM_THREADS=1
export RHOME=/project2/mstephens/software/R-3.5.3-mkl
export PATH=$RHOME/bin:$PATH
export LD_LIBRARY_PATH=$RHOME/lib64/R/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH=$RHOME/lib64/R/lib:$LIBRARY_PATH


###Loop over replicates
for (( rep=1; rep<=${n_reps}; rep++ ))
do

  ###Set seed for R scripts
  ranseed=`expr 32 + ${rep}`
  
  ###Make some necessary directories
  mkdir -p "${inter_outdir}/rep${rep}/prior/sumstats"
  mkdir -p "${inter_outdir}/rep${rep}/prior/matrices"
  mkdir -p "${inter_outdir}/rep${rep}/misc"

  cd ../dsc
  
  ################
  ####Step 1
  ################

  ###Sample training and test set indeces
  Rscript ../code/sample_testset_ids_gtex.R --output "${inter_outdir}/rep${rep}/misc/testset_ids.rds" \
                                            --n_testset "${n_testind}" \
                                            --seed "${ranseed}"
                                           

  ################
  ####Step 2
  ################

  ###Run dsc to simulate genes and process data for data-driven covariance estimation                                       
  "${dsc_script}" -o "${dsc_outdir}/rep${rep}" --target sim_proc --n_dataset "${total_genes}" -c 4 --host midway2.yml \
                  --testset_ids "\"${inter_outdir}/rep${rep}/misc/testset_ids.rds\"" \
                  --data_driven_mats_file "\"${inter_outdir}/rep${rep}/prior/matrices/${dsc_outname}.EZ.FL_PC3.rds\""
  

  ################
  ####Step 3
  ################

  ###Prepare summary stats for data-driven covariance estimation
  Rscript ../code/prepare_sumstats_for_ED_prior.R --input_path "${dsc_outdir}/rep${rep}/univ_sumstats" \
                                                  --output "${inter_outdir}/rep${rep}/prior/sumstats/${dsc_outname}_effects_for_ED_prior.rds" \
                                                  --n_random 4 \
                                                  --seed "${ranseed}"

  ###Run GTEx pipeline to compute the prior covariance matrices
  sos run /project2/mstephens/fmorgante/gtexresults/workflows/mashr_flashr_workflow.ipynb prior \
          -c ../scripts/midway2_sos.yml \
          -q midway2 \
          --cwd "${inter_outdir}/rep${rep}/prior/matrices" \
          --data "${inter_outdir}/rep${rep}/prior/sumstats/${dsc_outname}_effects_for_ED_prior.rds" \
          --output-prefix "${dsc_outname}"
          

  ################
  ####Step 4
  ################

  ###Run dsc to fit the models and do prediction
  "${dsc_script}" -o "${dsc_outdir}/rep${rep}" --target fit_pred_score -s existing --n_dataset "${pred_genes}" -c 4 --host midway2.yml \
                  --testset_ids "\"${inter_outdir}/rep${rep}/misc/testset_ids.rds\"" \
                  --data_driven_mats_file "\"${inter_outdir}/rep${rep}/prior/matrices/${dsc_outname}.EZ.FL_PC3.rds\""


  echo "Finished rep ${rep}"

done

###Deactivate the environment
source deactivate
