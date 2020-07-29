#!/usr/bin/env dsc

## A DSC for evaluating prediction accuracy of
## mr.mash in different scenarios.
DSC:
  R_libs: mr.mash.alpha, glmnet
  lib_path: functions
  exec_path: modules
  replicate: 2
  define:
    simulate: indepX_indepV_indepB_allr_norm, corrX_indepV_indepB, highcorrX_indepV_indepB, 
              indepX_indepV_sharedB, corrX_indepV_sharedB, highcorrX_indepV_sharedB
    fit:      mr_mash_em_singletons_no_datadriven, mr_mash_em_singletons_no_datadriven_drop_w0,
              mlasso, mridge, menet
    predict:  predict_linear
    score:    r2, scaled_mse, bias
  run: simulate * fit * predict * score

## Simulate modules
#Independent predictors, independent residuals, independent effects from a single normal,
#all resposens are causal
indepX_indepV_indepB_allr_norm: simulate_data_mod.R
  n:        900
  p:        50
  p_causal: 5
  r:        5
  r_causal: 5
  pve:      0.15
  B_cor:    0
  B_scale:  1
  w:        1
  X_cor:    0
  X_scale:  1
  V_cor:    0
  prop_testset: 0.2
  $Xtrain: out$Xtrain
  $Ytrain: out$Ytrain
  $Xtest:  out$Xtest
  $Ytest:  out$Ytest
  $B_true: out$B_true
 
#Correlated predictors, independent residuals, independent effects from a single normal,
#all resposens are causal
corrX_indepV_indepB_allr_norm(indepX_indepV_indepB_allr_norm):
  X_cor:   0.5
  
#Higly correlated predictors, independent residuals, independent effects from a single normal,
#all resposens are causal
highcorrX_indepV_indepB_allr_norm(indepX_indepV_indepB_allr_norm):
  X_cor:   0.8
  
#Independent predictors, independent residuals, shared effects from a single normal,
#all resposens are causal
indepX_indepV_sharedB_allr_norm(indepX_indepV_indepB_allr_norm):
  B_cor:   1

#Correlated predictors, independent residuals, shared effects from a single normal,
#all resposens are causal
corrX_indepV_sharedB_allr_norm(indepX_indepV_indepB_allr_norm):
  B_cor:   1
  X_cor:   0.5

#Higly correlated predictors, independent residuals, shared effects from a single normal,
#all resposens are causal
highcorrX_indepV_sharedB_allr_norm(indepX_indepV_indepB_allr_norm):
  B_cor:   1
  X_cor:   0.8


## Fit modules
#EM w0 updates, no drop w0, standardize X, update V (constrained diagonal),
#singletons, no data-driven matrices
mr_mash_em_singletons_no_datadriven: fit_mr_mash_mod.R
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
  singletons:             TRUE
  hetgrid:                (0, 0.5, 1)
  data_driven_mats:       FALSE
  subset_thresh:          0.05
  nthreads:               1
  $fit_obj:               out$fit
  $B_est:                 out$B_est
  $intercept_est:         out$intercept_est
  $time:                  out$elapsed_time
  
#EM w0 updates, drop w0 < 1e-8, standardize X, update V (constrained diagonal),
#singletons, no data-driven matrices
mr_mash_em_singletons_no_datadriven_drop_w0(mr_mash_em_singletons_no_datadriven):
  w0_threshold: 1e-8
  
#Multivariate LASSO  
mlasso: fit_mglmnet_mod.R
  X:                    $Xtrain
  Y:                    $Ytrain
  alpha:                1
  standardize:          TRUE
  $fit_obj:             out$fit
  $B_est:               out$B_est
  $intercept_est:       out$intercept_est
  $time:                out$elapsed_time

#Multivariate ridge  
mridge(mlasso):
  alpha:                0

#Multivariate enet  
menet(mlasso):
  alpha:                0.5


## Predict module
predict_linear: predict_mod.R
  B:         $B_est
  intercept: $intercept_est
  X:         $Xtest
  $Yhattest: Yhattest

## Score modules
r2: r2_mod.R
  Y:    $Ytest
  Yhat: $Yhattest 
  $err: err

scaled_mse: scaled_mse_mod.R
  Y:    $Ytest
  Yhat: $Yhattest 
  $err: err
  
bias: bias_mod.R
  Y:    $Ytest
  Yhat: $Yhattest 
  $err: err
  
