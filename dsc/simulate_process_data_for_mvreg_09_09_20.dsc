#!/usr/bin/env dsc

## A DSC for evaluating prediction accuracy of
## mr.mash in different scenarios.
DSC:
  R_libs: mr.mash.alpha
  lib_path: functions
  exec_path: modules
  replicate: 50
  define:
    simulate: indepX_indepV_sharedB_allr_norm
    process: univ_sumstats
  run: simulate * process

## Simulate modules
#Independent predictors, independent residuals, independent effects from a single normal,
#all resposens are causal
indepX_indepV_indepB_allr_norm: simulate_data_mod.R
  n:        90
  p:        500
  p_causal: 10
  r:        5
  r_causal: raw(list(1:5))
  pve:      0.5
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
  X_cor:    0.5
  
#Higly correlated predictors, independent residuals, independent effects from a single normal,
#all resposens are causal
highcorrX_indepV_indepB_allr_norm(indepX_indepV_indepB_allr_norm):
  X_cor:    0.8
  
#Independent predictors, independent residuals, shared effects from a single normal,
#all resposens are causal
indepX_indepV_sharedB_allr_norm(indepX_indepV_indepB_allr_norm):
  B_cor:    1

#Correlated predictors, independent residuals, shared effects from a single normal,
#all resposens are causal
corrX_indepV_sharedB_allr_norm(indepX_indepV_indepB_allr_norm):
  B_cor:    1
  X_cor:    0.5

#Higly correlated predictors, independent residuals, shared effects from a single normal,
#all resposens are causal
highcorrX_indepV_sharedB_allr_norm(indepX_indepV_indepB_allr_norm):
  B_cor:    1
  X_cor:    0.8
  
#Independent predictors, independent residuals, independent effects from a 2-component mixture
#of normals, all resposens are causal with a 2-block structure
indepX_indepV_indepB_2blocksr_norm(indepX_indepV_indepB_allr_norm):
  r_causal: raw(list(1:10,11:50))
  B_scale:  (1,1)
  B_cor:    (1,0.5)
  w:        (0.5,0.5)


#Compute univariate summary statistics
univ_sumstats: get_univ_sumstats_mod.R
  X:            $Xtrain
  Y:            $Ytrain
  zscores:      FALSE
  standardize:  TRUE
  nthreads:     1
  $sumstats:    out
  



