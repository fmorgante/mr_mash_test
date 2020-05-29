#!/usr/bin/env dsc

## A DSC for evaluating prediction accuracy of
## mr.mash in different scenarios.
DSC:
  R_libs: mr.mash.alpha, mr.ash.alpha, glmnet, gflasso
  lib_path: functions
  exec_path: modules
  replicate: 20
  define:
    simulate: indepX_indepV_indepB, corrX_indepV_indepB, highcorrX_indepV_indepB, 
              indepX_indepV_1respB, corrX_indepV_1respB, highcorrX_indepV_1respB,
              indepX_indepV_sharedB, corrX_indepV_sharedB, highcorrX_indepV_sharedB
    fit:      mr_mash_consec_em, mr_mash_consec_em_init_mlasso, 
              mr_mash_consec_em_drop_w0, mr_mash_consec_em_drop_w0_init_mlasso,
              mlasso, mridge, menet, gflasso
              #mr_mash_consec_mixsqp, mr_mash_declogBF_mixsqp, mr_mash_consec_em_init_indep, 
              #mr_mash_consec_em_daarem_init_indep, mr_mash_consec_mixsqp_init_indep,
              #mr_mash_consec_mixsqp_init_shared, mr_mash_consec_mixsqp_init_2pass,
              #mr_mash_consec_mixsqp_init_trueB, mr_mash_consec_mixsqp_init_mlasso
              #mr_mash_declogBF_em, mr_mash_declogBF_em_daarem, 
              #mr_mash_consec_em_init_shared, mr_mash_consec_em_daarem_init_shared, 
              #mr_mash_consec_em_init_2pass, mr_mash_consec_em_daarem_init_2pass, 
              #mr_mash_consec_em_init_trueB, mr_mash_consec_em_daarem_init_trueB,
              #, mr_mash_consec_em_daarem_init_mlasso,

    predict:  predict_linear
    score:    r2, mse, bias
  run: simulate * fit * predict * score

## Simulate modules
#Independent predictors, lowly correlated residuals, independent effects
indepX_indepV_indepB: simulate_data_mod.R
  n:        600
  p:        1000
  p_causal: 50
  r:        10
  r_causal: 10
  pve:      0.5
  B_cor:    0
  B_scale:  0.8
  X_cor:    0
  X_scale:  0.8
  V_cor:    0
  prop_testset: 0.2
  $Xtrain: out$Xtrain
  $Ytrain: out$Ytrain
  $Xtest:  out$Xtest
  $Ytest:  out$Ytest
  $B_true: out$B_true
 
#Correlated predictors, lowly correlated residuals, independent effects
corrX_indepV_indepB(indepX_indepV_indepB):
  B_cor:   0
  B_scale: 0.8
  X_cor:   0.5
  X_scale: 0.8
  V_cor:   0
  
#Higly correlated predictors, lowly correlated residuals, independent effects
highcorrX_indepV_indepB(indepX_indepV_indepB):
  B_cor:   0
  B_scale: 0.8
  X_cor:   0.8
  X_scale: 0.8
  V_cor:   0
  
#Independent predictors, lowly correlated residuals, shared effects
indepX_indepV_sharedB(indepX_indepV_indepB):
  B_cor:   1
  B_scale: 0.8
  X_cor:   0
  X_scale: 0.8
  V_cor:   0

#Correlated predictors, lowly correlated residuals, shared effects
corrX_indepV_sharedB(indepX_indepV_indepB):
  B_cor:   1
  B_scale: 0.8
  X_cor:   0.5
  X_scale: 0.8
  V_cor:   0

#Higly correlated predictors, lowly correlated residuals, shared effects
highcorrX_indepV_sharedB(indepX_indepV_indepB):
  B_cor:   1
  B_scale: 0.8
  X_cor:   0.8
  X_scale: 0.8
  V_cor:   0
  
#Independent predictors, lowly correlated residuals, effects in only 1 resp
indepX_indepV_1respB(indepX_indepV_indepB):
  r_causal: 1
  B_cor:   0
  B_scale: 0.8
  X_cor:   0
  X_scale: 0.8
  V_cor:   0

#Correlated predictors, lowly correlated residuals, effects in only 1 resp
corrX_indepV_1respB(indepX_indepV_indepB):
  r_causal: 1
  B_cor:   0
  B_scale: 0.8
  X_cor:   0.5
  X_scale: 0.8
  V_cor:   0

#Higly correlated predictors, lowly correlated residuals, effects in only 1 resp
highcorrX_indepV_1respB(indepX_indepV_indepB):
  r_causal: 1
  B_cor:   0
  B_scale: 0.8
  X_cor:   0.8
  X_scale: 0.8
  V_cor:   0


## Fit modules
#EM w0 updates, consecutive coordinate ascent updates
mr_mash_consec_em: fit_mr_mash_mod.R
  X:                    $Xtrain
  Y:                    $Ytrain
  update_w0:            TRUE
  update_w0_method:     "EM"
  w0_threshold:         0
  standardize:          TRUE
  update_V:             TRUE
  ca_update_order:      "consecutive"
  init_method:          "default"
  convergence_criterion: "ELBO"
  tol:                  1e-2
  B_true:               $B_true
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
  
#EM w0 updates, consecutive coordinate ascent updates, mu1 initilized by mlasso
mr_mash_consec_em_init_mlasso(mr_mash_consec_em):
  init_method: "mlasso"

#mixsqp w0 updates, consecutive coordinate ascent updates, mu1 initilized by mlasso
mr_mash_consec_mixsqp_init_mlasso(mr_mash_consec_em):
  update_w0_method: "mixsqp"
  init_method:    "mlasso"
  
#EM w0 updates, consecutive coordinate ascent updates, mu1 initilized by mlasso, daarem
mr_mash_consec_em_daarem_init_mlasso(mr_mash_consec_em):
  init_method: "mlasso"
  daarem: TRUE
  
#EM updates, consecutive coordinate ascent updates, drop w0 < 1e-8
mr_mash_consec_em_drop_w0(mr_mash_consec_em):
  w0_threshold: 1e-8
  
#EM w0 updates, consecutive coordinate ascent updates, drop w0 < 1e-8, mu1 initilized by mlasso
mr_mash_consec_em_drop_w0_init_mlasso(mr_mash_consec_em):
  init_method: "mlasso"
  w0_threshold: 1e-8

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
mridge(mlasso):
  alpha:                0

#Multivariate enet  
menet(mlasso):
  alpha:                0.5
  
#Multivariate Graph-Guided fused LASSO  
gflasso: fit_gflasso_mod.R
  X:                    $Xtrain
  Y:                    $Ytrain
  nCores:               1
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
  
