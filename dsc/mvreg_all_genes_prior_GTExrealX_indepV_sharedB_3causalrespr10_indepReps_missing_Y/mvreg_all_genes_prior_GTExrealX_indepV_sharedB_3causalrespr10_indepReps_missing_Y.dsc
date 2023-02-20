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
    simulate: indepV_sharedB_subsetcausalr
    process: univ_sumstats
    mr_mash_em_can_mlasso_out: mlasso_init * mr_mash_em_can_mlasso
    mr_mash_em_data_mlasso_out: mlasso_init * mr_mash_em_data_mlasso
    mr_mash_em_dataAndcan_mlasso_out: mlasso_init * mr_mash_em_dataAndcan_mlasso
    mr_mash_em_dataAndcan_dropcomp_mlasso_out: mlasso_init * mr_mash_em_dataAndcan_dropcomp_mlasso
    mr_mash_em_can_enet_out: enet_init * mr_mash_em_can_enet
    mr_mash_em_data_enet_out: enet_init * mr_mash_em_data_enet
    mr_mash_em_data_mean_impute_enet_out: enet_init * mr_mash_em_data_mean_impute_enet
    enet_out: enet_init * enet
    fit: mr_mash_em_dataAndcan_dropcomp_mlasso_out, mr_mash_em_dataAndcan_mlasso_out, 
         mr_mash_em_can_mlasso_out, mr_mash_em_can_enet_out, mr_mash_em_data_mlasso_out, 
         mr_mash_em_data_enet_out, mr_mash_em_data_mean_impute_enet_out, mtlasso, enet_out
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
indepV_sharedB_subsetcausalr: simulate_data_all_genes_prior_gtex_missing_Y_mod.R
  X:              $X
  p_causal:       "1:10"
  r:              10
  r_causal:       raw(list(1:3))
  B_scale:        1
  B_cor:          1
  w:              1
  pve:            0.2
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
  mu1_init:               NULL
  $fit_obj:               out$fit
  $B_est:                 out$B_est
  $intercept_est:         out$intercept_est
  $time:                  out$elapsed_time
  
#EM w0 updates, standardize X, update V (constrained diagonal),
#data-driven matrices, mlasso initialization
mr_mash_em_can_mlasso(mr_mash_em_can):
  mu1_init:               $B_est_init

#EM w0 updates, standardize X, update V (constrained diagonal),
#data-driven matrices, mlasso initialization
mr_mash_em_data_mlasso(mr_mash_em_can_mlasso):
  canonical_mats:         FALSE
  data_driven_mats:       ${data_driven_mats_file}

#EM w0 updates, standardize X, update V (constrained diagonal),
#canonical and data-driven matrices, mlasso initialization
mr_mash_em_dataAndcan_mlasso(mr_mash_em_data_mlasso):
  canonical_mats:         TRUE
  
#EM w0 updates, standardize X, update V (constrained diagonal),
#canonical and data-driven matrices, mlasso initialization
mr_mash_em_dataAndcan_dropcomp_mlasso(mr_mash_em_dataAndcan_mlasso):
  w0_threshold:           1e-08

#EM w0 updates, standardize X, update V (constrained diagonal),
#canonical matrices, enet initialization
mr_mash_em_can_enet(mr_mash_em_can):
  mu1_init:               $B_est_init

#EM w0 updates, standardize X, update V (constrained diagonal),
#data-driven matrices, enet initialization
mr_mash_em_data_enet(mr_mash_em_can):
  canonical_mats:         FALSE
  data_driven_mats:       ${data_driven_mats_file}
  mu1_init:               $B_est_init

#EM w0 updates, no drop w0, standardize X, update V (constrained diagonal),
#canonical matrices, mean impute Y
mr_mash_em_can_mean_impute: fit_mr_mash_all_genes_prior_mean_impute_mod.R
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
  mu1_init:               NULL
  $fit_obj:               out$fit
  $B_est:                 out$B_est
  $intercept_est:         out$intercept_est
  $time:                  out$elapsed_time

#EM w0 updates, standardize X, update V (constrained diagonal),
#data-driven matrices, mean imputed Y enet initialization
mr_mash_em_data_mean_impute_enet(mr_mash_em_can_mean_impute):
  canonical_mats:         FALSE
  data_driven_mats:       ${data_driven_mats_file}
  mu1_init:               $B_est_init

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
enet_init: fit_glmnet_missing_Y_mod.R
  X:                    $Xtrain
  Y:                    $Ytrain
  alpha:                0.5
  standardize:          TRUE
  nthreads:             1
  $B_est_init:          out$B_est
  $intercept_est_init:  out$intercept_est
  $time_init:           out$elapsed_time
  
enet: R(B_est <- B_in; intercept_est <- intercept_in; time <- time_in)
  B_in:                 $B_est_init
  intercept_in:         $intercept_est_init
  time_in:              $time_init
  $B_est:               B_est
  $intercept_est:       intercept_est
  $time:                time
  

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
  
