#!/usr/bin/env dsc

## A DSC for evaluating prediction accuracy of
## mr.mash in different scenarios.
DSC:
  R_libs: mr.mash.alpha, mr.ash.alpha
  lib_path: functions
  exec_path: modules
  global:
    n: 600
    r: 10
    p: 1000
    causal: 50 
    pve: 0.5
  replicate: 50
  define:
    simulate: indepX_lowcorrV_indepB, corrX_lowcorrV_indepB, 
              indepX_lowcorrV_sharedB, corrX_lowcorrV_sharedB             
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
  n: ${n}
  p: ${p}
  p_causal: ${causal}
  r: ${r}
  pve: ${pve}
  Sigma_cor_offdiag: 0
  Sigma_scale: 0.8
  Gamma_cor_offdiag: 0
  Gamma_scale: 0.8
  V_cor_offdiag: 0.15
  V_offdiag_scale: 1
  prop_testset: 0.2
  $Xtrain: out$Xtrain
  $Ytrain: out$Ytrain
  $Xtest: out$Xtest
  $Ytest: out$Ytest
  
#Correlated predictors, lowly correlated residuals, independent effects
corrX_lowcorrV_indepB: simulate_data_mod.R
  n: ${n}
  p: ${p}
  p_causal: ${causal}
  r: ${r}
  pve: ${pve}
  Sigma_cor_offdiag: 0
  Sigma_scale: 0.8
  Gamma_cor_offdiag: 0.5
  Gamma_scale: 0.8
  V_cor_offdiag: 0.15
  V_offdiag_scale: 1
  prop_testset: 0.2
  $Xtrain: out$Xtrain
  $Ytrain: out$Ytrain
  $Xtest: out$Xtest
  $Ytest: out$Ytest
  
#Independent predictors, lowly correlated residuals, shared effects
indepX_lowcorrV_sharedB: simulate_data_mod.R
  n: ${n}
  p: ${p}
  p_causal: ${causal}
  r: ${r}
  pve: ${pve}
  Sigma_cor_offdiag: 1
  Sigma_scale: 0.8
  Gamma_cor_offdiag: 0
  Gamma_scale: 0.8
  V_cor_offdiag: 0.15
  V_offdiag_scale: 1
  prop_testset: 0.2
  $Xtrain: out$Xtrain
  $Ytrain: out$Ytrain
  $Xtest: out$Xtest
  $Ytest: out$Ytest

#Correlated predictors, lowly correlated residuals, shared effects
corrX_lowcorrV_sharedB: simulate_data_mod.R
  n: ${n}
  p: ${p}
  p_causal: ${causal}
  r: ${r}
  pve: ${pve}
  Sigma_cor_offdiag: 1
  Sigma_scale: 0.8
  Gamma_cor_offdiag: 0.5
  Gamma_scale: 0.8
  V_cor_offdiag: 0.15
  V_offdiag_scale: 1
  prop_testset: 0.2
  $Xtrain: out$Xtrain
  $Ytrain: out$Ytrain
  $Xtest: out$Xtest
  $Ytest: out$Ytest


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
  $fit_obj: out$fit
  $time: out$elapsed_time
  
#mixsqp w0 updates, consecutive coordinate ascent updates
mr_mash_consec_mixsqp: fit_mr_mash_mod.R
  X: $Xtrain
  Y: $Ytrain
  update_w0: TRUE
  update_w0_method: "mixsqp"
  standardize: TRUE
  update_V: TRUE
  ca_update_order: "consecutive"
  mr_ash_method: NULL
  $fit_obj: out$fit
  $time: out$elapsed_time

#EM w0 updates, decreasing logBF coordinate ascent updates
mr_mash_declogBF_em: fit_mr_mash_mod.R
  X: $Xtrain
  Y: $Ytrain
  update_w0: TRUE
  update_w0_method: "EM"
  standardize: TRUE
  update_V: TRUE
  ca_update_order: "decreasing_logBF"
  mr_ash_method: NULL
  $fit_obj: out$fit
  $time: out$elapsed_time
  
#mixsqp w0 updates, decreasing logBF coordinate ascent updates
mr_mash_declogBF_mixsqp: fit_mr_mash_mod.R
  X: $Xtrain
  Y: $Ytrain
  update_w0: TRUE
  update_w0_method: "mixsqp"
  standardize: TRUE
  update_V: TRUE
  ca_update_order: "decreasing_logBF"
  mr_ash_method: NULL
  $fit_obj: out$fit
  $time: out$elapsed_time
  
#EM w0 updates, consecutive coordinate ascent updates, mu1 initilized by mr.ash
#run on each response
mr_mash_consec_em_init_indep: fit_mr_mash_mod.R
  X: $Xtrain
  Y: $Ytrain
  update_w0: TRUE
  update_w0_method: "EM"
  standardize: TRUE
  update_V: TRUE
  ca_update_order: "consecutive"
  mr_ash_method: "independent"
  $fit_obj: out$fit
  $time: out$elapsed_time

#mixsqp w0 updates, consecutive coordinate ascent updates, mu1 initilized by mr.ash
#run on each response
mr_mash_consec_mixsqp_init_indep: fit_mr_mash_mod.R
  X: $Xtrain
  Y: $Ytrain
  update_w0: TRUE
  update_w0_method: "mixsqp"
  standardize: TRUE
  update_V: TRUE
  ca_update_order: "consecutive"
  mr_ash_method: "independent"
  $fit_obj: out$fit
  $time: out$elapsed_time
  
#EM w0 updates, consecutive coordinate ascent updates, mu1 initilized by mr.ash
#run on stacked responses
mr_mash_consec_em_init_shared: fit_mr_mash_mod.R
  X: $Xtrain
  Y: $Ytrain
  update_w0: TRUE
  update_w0_method: "EM"
  standardize: TRUE
  update_V: TRUE
  ca_update_order: "consecutive"
  mr_ash_method: "shared"
  $fit_obj: out$fit
  $time: out$elapsed_time

#mixsqp w0 updates, consecutive coordinate ascent updates, mu1 initilized by mr.ash
#run on stacked responses
mr_mash_consec_mixsqp_init_shared: fit_mr_mash_mod.R
  X: $Xtrain
  Y: $Ytrain
  update_w0: TRUE
  update_w0_method: "mixsqp"
  standardize: TRUE
  update_V: TRUE
  ca_update_order: "consecutive"
  mr_ash_method: "shared"
  $fit_obj: out$fit
  $time: out$elapsed_time
  

## Predict module
predict_linear: predict_mod.R
  fit: $fit_obj
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
  
