#!/usr/bin/env dsc

## A DSC for evaluating prediction accuracy of
## mr.mash in different scenarios.
DSC:
  R_libs: mr.mash.alpha, mr.ash.alpha
  lib_path: functions
  exec_path: modules
  replicate: 20
  define:
    simulate: indepX_lowcorrV_indepB, corrX_lowcorrV_indepB, highcorrX_lowcorrV_indepB, 
              indepX_lowcorrV_sharedB, corrX_lowcorrV_sharedB, highcorrX_lowcorrV_sharedB             
    fit:      mr_mash_consec_em, mr_mash_consec_em_daarem,  
              mr_mash_declogBF_em, mr_mash_declogBF_em_daarem, 
              mr_mash_consec_em_init_shared, mr_mash_consec_em_daarem_init_shared, 
              mr_mash_consec_em_init_2pass, mr_mash_consec_em_daarem_init_2pass, 
              mr_mash_consec_em_init_trueB, mr_mash_consec_em_daarem_init_trueB,
              mlasso, mridge, menet
              #mr_mash_consec_mixsqp, mr_mash_declogBF_mixsqp, mr_mash_consec_em_init_indep, 
              #mr_mash_consec_em_daarem_init_indep, mr_mash_consec_mixsqp_init_indep,
              #mr_mash_consec_mixsqp_init_shared, mr_mash_consec_mixsqp_init_2pass,
              #mr_mash_consec_mixsqp_init_trueB
    predict:  predict_linear
    score:    r2, mse, bias
  run: simulate * fit * predict * score

## Simulate modules
#Independent predictors, lowly correlated residuals, independent effects
indepX_lowcorrV_indepB: simulate_data_mod.R
  n:        600
  p:        1000
  p_causal: 50, 500
  r:        10
  pve:      0.5
  B_cor:    0
  B_scale:  0.8
  X_cor:    0
  X_scale:  0.8
  V_cor:    0.15
  prop_testset: 0.2
  $Xtrain: out$Xtrain
  $Ytrain: out$Ytrain
  $Xtest:  out$Xtest
  $Ytest:  out$Ytest
  $B_true: out$B_true
 
#Correlated predictors, lowly correlated residuals, independent effects
corrX_lowcorrV_indepB(indepX_lowcorrV_indepB):
  B_cor:   0
  B_scale: 0.8
  X_cor:   0.5
  X_scale: 0.8
  V_cor:   0.15
  
#Higly correlated predictors, lowly correlated residuals, independent effects
highcorrX_lowcorrV_indepB(indepX_lowcorrV_indepB):
  B_cor:   0
  B_scale: 0.8
  X_cor:   0.8
  X_scale: 0.8
  V_cor:   0.15
  
#Independent predictors, lowly correlated residuals, shared effects
indepX_lowcorrV_sharedB(indepX_lowcorrV_indepB):
  B_cor:   1
  B_scale: 0.8
  X_cor:   0
  X_scale: 0.8
  V_cor:   0.15

#Correlated predictors, lowly correlated residuals, shared effects
corrX_lowcorrV_sharedB(indepX_lowcorrV_indepB):
  B_cor:   1
  B_scale: 0.8
  X_cor:   0.5
  X_scale: 0.8
  V_cor:   0.15

#Higly correlated predictors, lowly correlated residuals, shared effects
highcorrX_lowcorrV_sharedB(indepX_lowcorrV_indepB):
  B_cor:   1
  B_scale: 0.8
  X_cor:   0.8
  X_scale: 0.8
  V_cor:   0.15

## Fit modules
#EM w0 updates, consecutive coordinate ascent updates
mr_mash_consec_em: fit_mr_mash_mod.R
  X:                    $Xtrain
  Y:                    $Ytrain
  update_w0:            TRUE
  update_w0_method:     "EM"
  standardize:          TRUE
  update_V:             TRUE
  ca_update_order:      "consecutive"
  init_method:          "default"
  B_true:               $B_true
  select_w0_threshold:  0
  daarem:               FALSE
  $fit_obj:             out$fit
  $B_est:               out$B_est
  $intercept_est:       out$intercept_est
  $time:                out$elapsed_time
  
#mixsqp w0 updates, consecutive coordinate ascent updates
mr_mash_consec_mixsqp(mr_mash_consec_em):
  update_w0_method: "mixsqp"

#EM updates, consecutive coordinate ascent updates, daarem
mr_mash_consec_em_daarem(mr_mash_consec_em):
  daarem: TRUE

#EM w0 updates, decreasing logBF coordinate ascent updates
mr_mash_declogBF_em(mr_mash_consec_em):
  ca_update_order: "decreasing_logBF"
  
#mixsqp w0 updates, decreasing logBF coordinate ascent updates
mr_mash_declogBF_mixsqp(mr_mash_consec_em):
  update_w0_method: "mixsqp"
  ca_update_order:  "decreasing_logBF"
  
#EM updates, decreasing logBF coordinate ascent updates, daarem
mr_mash_declogBF_em_daarem(mr_mash_consec_em):
  ca_update_order:  "decreasing_logBF"
  daarem: TRUE
  
#EM w0 updates, consecutive coordinate ascent updates, mu1 initilized by mr.ash
#run on each response
mr_mash_consec_em_init_indep(mr_mash_consec_em):
  init_method: "independent"

#mixsqp w0 updates, consecutive coordinate ascent updates, mu1 initilized by mr.ash
#run on each response
mr_mash_consec_mixsqp_init_indep(mr_mash_consec_em):
  update_w0_method: "mixsqp"
  init_method:    "independent"
  
#EM w0 updates, consecutive coordinate ascent updates, mu1 initilized by mr.ash
#run on each response, daarem
mr_mash_consec_em_daarem_init_indep(mr_mash_consec_em):
  init_method:    "independent"
  daarem: TRUE
  
#EM w0 updates, consecutive coordinate ascent updates, mu1 initilized by mr.ash
#run on stacked responses
mr_mash_consec_em_init_shared(mr_mash_consec_em):
  init_method: "shared"

#mixsqp w0 updates, consecutive coordinate ascent updates, mu1 initilized by mr.ash
#run on stacked responses
mr_mash_consec_mixsqp_init_shared(mr_mash_consec_em):
  update_w0_method: "mixsqp"
  init_method:    "shared"
  
#EM w0 updates, consecutive coordinate ascent updates, mu1 initilized by mr.ash
#run on stacked responses, daarem
mr_mash_consec_em_daarem_init_shared(mr_mash_consec_em):
  init_method:    "shared"
  daarem: TRUE

#EM w0 updates, consecutive coordinate ascent updates, mu1 initilized by mr.ash
#two pass
mr_mash_consec_em_init_2pass(mr_mash_consec_em):
  init_method: "2pass"

#mixsqp w0 updates, consecutive coordinate ascent updates, mu1 initilized by mr.ash
#two pass
mr_mash_consec_mixsqp_init_2pass(mr_mash_consec_em):
  update_w0_method: "mixsqp"
  init_method:    "2pass"
  
#EM w0 updates, consecutive coordinate ascent updates, mu1 initilized by mr.ash
#two pass, daarem
mr_mash_consec_em_daarem_init_2pass(mr_mash_consec_em):
  init_method: "2pass"
  daarem: TRUE
  
#EM w0 updates, consecutive coordinate ascent updates, mu1 initilized by true B
mr_mash_consec_em_init_trueB(mr_mash_consec_em):
  init_method: "truth"

#mixsqp w0 updates, consecutive coordinate ascent updates, mu1 initilized by true B
mr_mash_consec_mixsqp_init_trueB(mr_mash_consec_em):
  update_w0_method: "mixsqp"
  init_method:    "truth"
  
#EM w0 updates, consecutive coordinate ascent updates, mu1 initilized by true B, daarem
mr_mash_consec_em_daarem_init_trueB(mr_mash_consec_em):
  init_method: "truth"
  daarem: TRUE

#Multivariate LASSO  
mlasso: fit_mglmnet_mod.R
  X:                    $Xtrain
  Y:                    $Ytrain
  alpha:                1
  $fit_obj:             out$fit
  $B_est:               out$B_est
  $intercept_est:       out$intercept_est
  $time:                out$elapsed_time

#Multivariate ridge  
mridge: fit_mglmnet_mod.R
  X:                    $Xtrain
  Y:                    $Ytrain
  alpha:                0
  $fit_obj:             out$fit
  $B_est:               out$B_est
  $intercept_est:       out$intercept_est
  $time:                out$elapsed_time

#Multivariate enet  
menet: fit_mglmnet_mod.R
  X:                    $Xtrain
  Y:                    $Ytrain
  alpha:                0.5
  $fit_obj:             out$fit
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
r2: r2_mod.R
  Y:    $Ytest
  Yhat: $Yhattest 
  $err: err

mse: mse_mod.R
  Y:    $Ytest
  Yhat: $Yhattest 
  $err: err
  
bias: bias_mod.R
  Y:    $Ytest
  Yhat: $Yhattest 
  $err: err
  
