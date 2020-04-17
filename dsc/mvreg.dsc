#!/usr/bin/env dsc

## A DSC for evaluating prediction accuracy of
## mr.mash in different scenarios.
DSC:
  R_libs: mr.mash.alpha, mr.ash.alpha
  lib_path: functions
  exec_path: modules
  replicate: 50
  define:
    simulate: indepX_lowcorrV_indepB, corrX_lowcorrV_indepB, highcorrX_lowcorrV_indepB, 
              indepX_lowcorrV_sharedB, corrX_lowcorrV_sharedB, highcorrX_lowcorrV_sharedB             
    fit:      mr_mash_consec_em, mr_mash_consec_mixsqp, 
              mr_mash_declogBF_em, mr_mash_declogBF_mixsqp,
              mr_mash_consec_em_init_indep, mr_mash_consec_mixsqp_init_indep,
              mr_mash_consec_em_init_shared, mr_mash_consec_mixsqp_init_shared
    predict:  predict_linear
    score:    r2, mse, bias
  run: simulate * fit * predict * score

## Simulate modules
#Independent predictors, lowly correlated residuals, independent effects
indepX_lowcorrV_indepB: simulate_data_mod.R
  n: 600
  p: 1000
  p_causal: 50
  r: 10
  intercepts: R{rep(1, 10)}
  pve: 0.5
  B_cor: 0
  B_scale: 0.8
  X_cor: 0
  X_scale: 0.8
  V_cor: 0.15
  prop_testset: 0.2
  $Xtrain: out$Xtrain
  $Ytrain: out$Ytrain
  $Xtest: out$Xtest
  $Ytest: out$Ytest
  
#Correlated predictors, lowly correlated residuals, independent effects
corrX_lowcorrV_indepB(indepX_lowcorrV_indepB):
  B_cor: 0
  B_scale: 0.8
  X_cor: 0.5
  X_scale: 0.8
  V_cor: 0.15
  
#Higly correlated predictors, lowly correlated residuals, independent effects
highcorrX_lowcorrV_indepB(indepX_lowcorrV_indepB):
  B_cor: 0
  B_scale: 0.8
  X_cor: 0.8
  X_scale: 0.8
  V_cor: 0.15
  
#Independent predictors, lowly correlated residuals, shared effects
indepX_lowcorrV_sharedB(indepX_lowcorrV_indepB):
  B_cor: 1
  B_scale: 0.8
  X_cor: 0
  X_scale: 0.8
  V_cor: 0.15

#Correlated predictors, lowly correlated residuals, shared effects
corrX_lowcorrV_sharedB(indepX_lowcorrV_indepB):
  B_cor: 1
  B_scale: 0.8
  X_cor: 0.5
  X_scale: 0.8
  V_cor: 0.15

#Higly correlated predictors, lowly correlated residuals, shared effects
highcorrX_lowcorrV_sharedB(indepX_lowcorrV_indepB):
  B_cor: 1
  B_scale: 0.8
  X_cor: 0.8
  X_scale: 0.8
  V_cor: 0.15


## Fit modules
#EM w0 updates, consecutive coordinate ascent updates
mr_mash_consec_em: fit_mr_mash_mod.R
  X: $Xtrain
  Y: $Ytrain
  update_w0: TRUE
  update_w0_method: "EM"
  standardize: TRUE
  update_V: TRUE
  ca_update_order: "consecutive"
  mr_ash_method: NULL
  scaling_grid: R{seq(0.1, 2.1, 0.2)}
  $fit_obj: out$fit
  $B_est: out$B_est
  $intercept_est: out$intercept_est
  $time: out$elapsed_time
  
#mixsqp w0 updates, consecutive coordinate ascent updates
mr_mash_consec_mixsqp(mr_mash_consec_em):
  update_w0_method: "mixsqp"

#EM w0 updates, decreasing logBF coordinate ascent updates
mr_mash_declogBF_em(mr_mash_consec_em):
  ca_update_order: "decreasing_logBF"
  
#mixsqp w0 updates, decreasing logBF coordinate ascent updates
mr_mash_declogBF_mixsqp(mr_mash_consec_em):
  update_w0_method: "mixsqp"
  ca_update_order: "decreasing_logBF"
  
#EM w0 updates, consecutive coordinate ascent updates, mu1 initilized by mr.ash
#run on each response
mr_mash_consec_em_init_indep(mr_mash_consec_em):
  mr_ash_method: "independent"

#mixsqp w0 updates, consecutive coordinate ascent updates, mu1 initilized by mr.ash
#run on each response
mr_mash_consec_mixsqp_init_indep(mr_mash_consec_em):
  update_w0_method: "mixsqp"
  mr_ash_method: "independent"
  
#EM w0 updates, consecutive coordinate ascent updates, mu1 initilized by mr.ash
#run on stacked responses
mr_mash_consec_em_init_shared(mr_mash_consec_em):
  mr_ash_method: "shared"

#mixsqp w0 updates, consecutive coordinate ascent updates, mu1 initilized by mr.ash
#run on stacked responses
mr_mash_consec_mixsqp_init_shared(mr_mash_consec_em):
  update_w0_method: "mixsqp"
  mr_ash_method: "shared"
  

## Predict module
predict_linear: predict_mod.R
  B: $B_est
  intercept: $intercept_est
  X: $Xtest
  $Yhattest: Yhattest
  

## Score modules
r2: r2_mod.R
  Y: $Ytest
  Yhat: $Yhattest 
  $err: err

mse: mse_mod.R
  Y: $Ytest
  Yhat: $Yhattest 
  $err: err
  
bias: bias_mod.R
  Y: $Ytest
  Yhat: $Yhattest 
  $err: err
  
