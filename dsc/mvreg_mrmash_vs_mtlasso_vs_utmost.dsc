#!/usr/bin/env dsc

## A DSC for evaluating prediction accuracy of
## mr.mash in different scenarios.
DSC:
  R_libs: mr.mash.alpha, glmnet
  python_modules: numpy, mtlasso, pandas, sklearn.model_selection, time
  lib_path: functions
  exec_path: modules
  replicate: 1
  define:
    simulate: indepX_indepV_sharedB_allr_norm, indepX_indepV_indepB_allr_norm,
              indepX_indepV_B_1causalr, indepX_indepV_sharedB_2blocksr_diffPVE,
              indepX_indepV_sharedB_3causalrespr10
    fit:      mr_mash, mtlasso, utmost
    score:    mse
  run: simulate * fit * score

## Simulate modules
#Independent predictors, independent residuals, shared effects from a single normal,
#all resposens are causal
indepX_indepV_sharedB_allr_norm: simulate_data_mod.R
  n:        500
  p:        1000
  p_causal: "1:10"
  r:        10
  r_causal: raw(list(1:10))
  pve:      0.2
  B_cor:    1
  B_scale:  1
  w:        1
  X_cor:    0
  X_scale:  1
  V_cor:    0
  prop_testset: 0
  $X: out$X
  $Y: out$Y
  $B_true: out$B_true
 

#Independent predictors, independent residuals, independent effects from a single normal,
#all resposens are causal
indepX_indepV_indepB_allr_norm(indepX_indepV_sharedB_allr_norm):
  B_cor:    0

#Independent predictors, independent residuals, effects only in tissue 1
indepX_indepV_B_1causalr(indepX_indepV_sharedB_allr_norm):
  r_causal: raw(list(1))
  B_cor:    0

#Independent predictors, independent residuals, 2 blocks of effects, with
#shared effects within block, no sharing between blocks, and different PVE
#for each block
indepX_indepV_sharedB_2blocksr_diffPVE(indepX_indepV_sharedB_allr_norm):
  r_causal: raw(list(1:3,4:10))
  B_scale:  (0.8,1)
  B_cor:    (0.9,0.75)
  w:        (0.5,0.5)
  pve:      raw(c(rep(0.2, 3), rep(0.05, 7)))

#Independent predictors, independent residuals, shared effects only in tissue 1-3
indepX_indepV_sharedB_3causalrespr10(indepX_indepV_sharedB_allr_norm):
  r_causal: raw(list(1:3))
  B_cor:    1
  
## Fit modules
#EM w0 updates, no drop w0, standardize X, update V (constrained diagonal),
#singletons, no data-driven matrices
mr_mash: fit_mr_mash_mod.R
  X:                      $X
  Y:                      $Y
  update_w0:              TRUE
  update_w0_method:       "EM"
  w0_threshold:           0
  standardize:            FALSE
  update_V:               TRUE
  update_V_method:        "diagonal"
  ca_update_order:        "consecutive"
  convergence_criterion:  "ELBO"
  tol:                    1e-2
  singletons:             TRUE
  hetgrid:                (0,0.25,0.5,0.75,1)
  data_driven_mats:       FALSE
  subset_thresh:          0
  nthreads:               1
  $fit_obj:               out$fit
  $B_est:                 out$B_est
  $intercept_est:         out$intercept_est
  $time:                  out$elapsed_time
  
#Sparse multi-task lasso  
mtlasso: fit_mtlasso_mod.py
  X:                    $X
  Y:                    $Y
  standardize:          False
  nfolds:               5
  B_init:               None
  grid_limits:          None
  grid_length:          10
  $B_est:               B_est
  $intercept_est:       intercept_est
  $time:                elapsed_time

#UTMOST
utmost: fit_utmost_mod.R
  X:                    $X
  Y:                    $Y
  nfolds:               5
  $B_est:               out$B_est
  $time:                out$elapsed_time

## Score modules
mse: coeff_mse_mod.R
  B:                    $B_true
  Bhat:                 $B_est
  $err:                 err
