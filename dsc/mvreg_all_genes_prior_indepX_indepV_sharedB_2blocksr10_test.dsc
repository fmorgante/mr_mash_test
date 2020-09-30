#!/usr/bin/env dsc

## A DSC for evaluating prediction accuracy of
## mr.mash in different scenarios.
DSC:
  R_libs: mr.mash.alpha, glmnet, doMC
  python_modules: numpy, mtlasso, pandas, sklearn.model_selection, time
  lib_path: functions
  exec_path: modules, functions
  define:
    simulate: indepX_indepV_sharedB_2blocksr
    process: univ_sumstats
    mr_mash_em_can_mlasso: mlasso_init * mr_mash_em_can
    mr_mash_em_data_mlasso: mlasso_init * mr_mash_em_data
    mr_mash_em_dataAndcan_mlasso: mlasso_init * mr_mash_em_dataAndcan
    mlasso_out: mlasso_init * mlasso
    fit: mr_mash_em_dataAndcan_mlasso, mr_mash_em_can_mlasso, mr_mash_em_data_mlasso,
         mlasso_out, mridge, menet, mtlasso
    predict:  predict_linear
    score:    r2, scaled_mse, bias
  run: 
    sim_proc: simulate * process
    fit_pred_score: simulate * process * fit * predict * score

## Simulate modules
#Independent predictors, independent residuals, shared effects from a 2-component mixture
#of normals, all resposens are causal with a 2-block structure
indepX_indepV_sharedB_2blocksr: simulate_data_all_genes_prior_mod.R
  n:              90
  p:              500
  p_causal:       5
  r:              10
  r_causal:       raw(list(1:3,4:10))
  B_scale:        (0.8,1)
  B_cor:          (1,1)
  w:              (0.5,0.5)
  pve:            0.2
  X_cor:          0
  X_scale:        1
  V_cor:          0
  testset_index:  "../output/mvreg_all_genes_prior_indepX_indepV_sharedB_2blocksr10_test_inter/misc/testset_indeces.rds"
  $Xtrain: out$Xtrain
  $Ytrain: out$Ytrain
  $Xtest:  out$Xtest
  $Ytest:  out$Ytest


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
  data_driven_mats:       "/project2/mstephens/fmorgante/mr_mash_test/output/mvreg_all_genes_prior_indepX_indepV_sharedB_2blocksr10_test_inter/prior/matrices/mvreg_all_genes_prior_indepX_indepV_sharedB_2blocksr10_test.EZ.FL_PC3.rds"

#EM w0 updates, standardize X, update V (constrained diagonal),
#canonical and data-driven matrices
mr_mash_em_dataAndcan(mr_mash_em_can):
  data_driven_mats:       "/project2/mstephens/fmorgante/mr_mash_test/output/mvreg_all_genes_prior_indepX_indepV_sharedB_2blocksr10_test_inter/prior/matrices/mvreg_all_genes_prior_indepX_indepV_sharedB_2blocksr10_test.EZ.FL_PC3.rds"

#Multivariate LASSO  
mlasso_init: fit_mglmnet_mod.R
  X:                    $Xtrain
  Y:                    $Ytrain
  alpha:                1
  standardize:          TRUE
  nthreads:             1
  $fit_obj:             out$fit
  $B_est_init:          out$B_est
  $intercept_est_init:  out$intercept_est
  $time_init:           out$elapsed_time
  
mlasso: R(B_est <- B_in; intercept_est <- intercept_in; time <- time_in)
  B_in:                 $B_est_init
  intercept_in:         $intercept_est_init
  time_in:              $time_init
  $B_est:               B_est
  $intercept_est:       intercept_est
  $time:                time

#Multivariate ridge  
mridge: fit_mglmnet_mod.R
  X:                    $Xtrain
  Y:                    $Ytrain
  alpha:                0
  standardize:          TRUE
  nthreads:             1
  $fit_obj:             out$fit
  $B_est:               out$B_est
  $intercept_est:       out$intercept_est
  $time:                out$elapsed_time

#Multivariate enet  
menet(mridge):
  alpha:                0.5
  
#Sparse multi-task LASSO (aka UTMOST)  
mtlasso: fit_mtlasso.py + fit_mtlasso_mod.py
  X:                    $Xtrain
  Y:                    $Ytrain
  standardize:          False
  nfolds:               5
  $B_est:               B_est
  $intercept_est:       intercept_est
  $time:                elapsed_time


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
scaled_mse: scaled_rmse_mod.R
  Y:    $Ytest
  Yhat: $Yhattest 
  $err: err

#Slope of y~yhat  
bias: bias_mod.R
  Y:    $Ytest
  Yhat: $Yhattest 
  $err: err
  