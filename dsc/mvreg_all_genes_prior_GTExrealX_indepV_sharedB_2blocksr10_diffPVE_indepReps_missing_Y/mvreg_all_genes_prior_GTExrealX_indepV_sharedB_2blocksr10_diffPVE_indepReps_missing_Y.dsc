#!/usr/bin/env dsc

## A DSC for evaluating prediction accuracy of
## mr.mash in different scenarios.
DSC:
  R_libs: mr.mash.alpha, glmnet, doMC, matrixStats, mvtnorm, MBSP
  python_modules: numpy, mtlasso, pandas, sklearn.model_selection, time
  lib_path: ../../functions
  exec_path: ../../modules
  global:
    data_file: ../../../data/gtex-v8-manifest-full-X.txt
    n_dataset: 1
    randseed: 1
    testset_ids: "/path/to/file"
    data_driven_mats_file: "/path/to/file"
  define:
    data: extract_X
    simulate: indepV_sharedB_2blocksr_diffPVE
    process: univ_sumstats
    mr_mash_em_can_mlasso: mlasso_init * mr_mash_em_can
    mr_mash_em_data_mlasso: mlasso_init * mr_mash_em_data
    mr_mash_em_dataAndcan_mlasso: mlasso_init * mr_mash_em_dataAndcan
    mr_mash_em_dataAndcan_dropcomp_mlasso: mlasso_init * mr_mash_em_dataAndcan_dropcomp
    fit: mr_mash_em_dataAndcan_dropcomp_mlasso,mr_mash_em_dataAndcan_mlasso, 
         mr_mash_em_can_mlasso, mr_mash_em_data_mlasso, mtlasso, enet
    predict: predict_linear
    score: scaled_rmse
  run: 
    sim_proc: data * simulate * process
    fit_pred_score: data * simulate * process * fit * predict * score
    
    
##Data module
#Extract X for specified number of randomly selected genes
extract_X: ../../functions/utils.R + R(data = readRDS(paste0("../../", dataset));
            X = filter_X(data$X, missing_rate_cutoff, maf_cutoff, var_cutoff))
  dataset: Shell{get_seeded_random() { seed="$1" ; openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt </dev/zero 2>/dev/null ; } ; 
    shuf -n ${n_dataset} ${data_file} --random-source=<(get_seeded_random ${randseed}) | sed 's/^/"/;s/$/"/'}
  missing_rate_cutoff: 0.05
  maf_cutoff: 0.05
  var_cutoff: 0.05
  $X: X


## Simulate modules
#Independent residuals, independent effects from a normal
indepV_sharedB_2blocksr_diffPVE: simulate_data_all_genes_prior_gtex_missing_Y_mod.R
  X:              $X
  p_causal:       "1:10"
  r:              10
  r_causal:       raw(list(1:3,4:10))
  B_scale:        (0.8,1)
  B_cor:          (0.9,0.75)
  w:              (0.5,0.5)
  pve:            raw(c(rep(0.2, 3), rep(0.05, 7)))
  V_cor:          0
  testset_index:  ${testset_ids}
  prop_miss:      0.7
  $Xtrain: out$Xtrain
  $Ytrain: out$Ytrain
  $Xtest:  out$Xtest
  $Ytest:  out$Ytest
  $B_true: out$B_true


##Process modules
#Compute univariate summary statistics
univ_sumstats: get_univ_sumstats_mod.R
  X:            $Xtrain
  Y:            $Ytrain
  zscores:      FALSE
  standardize:  TRUE
  nthreads:     1
  $sumstats:    out
  

## Fit modules
#EM w0 updates, no drop w0, standardize X, update V (constrained diagonal),
#canonical matrices
mr_mash_em_can: fit_mr_mash_all_genes_prior_mod.R
  X:                      $Xtrain
  Y:                      $Ytrain
  update_w0:              TRUE
  update_w0_method:       "EM"
  w0_threshold:           0
  standardize:            TRUE
  update_V:               TRUE
  update_V_method:        "diagonal"
  ca_update_order:        "consecutive"
  convergence_criterion:  "ELBO"
  tol:                    1e-2
  canonical_mats:         TRUE
  singletons:             TRUE
  hetgrid:                (0, 0.25, 0.5, 0.75, 1)
  sumstats:               $sumstats
  data_driven_mats:       NULL
  nthreads:               1
  mu1_init:               $B_est_init
  $fit_obj:               out$fit
  $B_est:                 out$B_est
  $intercept_est:         out$intercept_est
  $time:                  out$elapsed_time

#EM w0 updates, standardize X, update V (constrained diagonal),
#data-driven matrices
mr_mash_em_data(mr_mash_em_can):
  canonical_mats:         FALSE
  data_driven_mats:       ${data_driven_mats_file}

#EM w0 updates, standardize X, update V (constrained diagonal),
#canonical and data-driven matrices
mr_mash_em_dataAndcan(mr_mash_em_data):
  canonical_mats:         TRUE
  
#EM w0 updates, standardize X, update V (constrained diagonal),
#canonical and data-driven matrices
mr_mash_em_dataAndcan_dropcomp(mr_mash_em_dataAndcan):
  w0_threshold:           1e-08
  
#Multivariate LASSO estimates  
mlasso_init: compute_coefficients_mlasso_missing_Y_mod.R
  X:                    $Xtrain
  Y:                    $Ytrain
  standardize:          TRUE
  nthreads:             1
  $B_est_init:          out$B
  
#Sparse multi-task LASSO (aka UTMOST)  
mtlasso: fit_mtlasso_missing_Y_mod.py
  X:                    $Xtrain
  Y:                    $Ytrain
  standardize:          True
  nfolds:               5
  B_init:               None
  grid_limits:          None
  grid_length:          10
  $B_est:               B_est
  $intercept_est:       intercept_est
  $time:                elapsed_time
  
#Univariate enet  
enet: fit_glmnet_missing_Y_mod.R
  X:                    $Xtrain
  Y:                    $Ytrain
  alpha:                0.5
  standardize:          TRUE
  nthreads:             1
  $B_est:               out$B_est
  $intercept_est:       out$intercept_est
  $time:                out$elapsed_time
  

## Predict module
predict_linear: predict_mod.R
  B:         $B_est
  intercept: $intercept_est
  X:         $Xtest
  $Yhattest: Yhattest


## Score modules
#r^2
r2: r2_mod.R
  Y:    $Ytest
  Yhat: $Yhattest 
  $err: err

#RMSE scaled by sd(y)
scaled_rmse: scaled_rmse_mod.R
  Y:    $Ytest
  Yhat: $Yhattest 
  $err: err

#Slope of y~yhat  
bias: bias_mod.R
  Y:    $Ytest
  Yhat: $Yhattest 
  $err: err
  