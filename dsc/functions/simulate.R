simulate_data <- function(n, p, p_causal, r, intercepts=rep(1, r),
                          pve, Sigma_cor_offdiag, Sigma_scale,
                          Gamma_cor_offdiag, Gamma_scale,
                          V_cor_offdiag, V_offdiag_scale, prop_testset){
  
  dat <- mr.mash.alpha:::simulate_mr_mash_data(n, p, p_causal, r,
                                               pve=pve, Sigma_cor_offdiag=Sigma_cor_offdiag, Sigma_scale=Sigma_scale,
                                               Gamma_cor_offdiag=Gamma_cor_offdiag, Gamma_scale=Gamma_scale,
                                               V_cor_offdiag=V_cor_offdiag, V_offdiag_scale=V_offdiag_scale)
  
  test_set <- sort(sample(x=c(1:n), size=round(n*prop_testset), replace=FALSE))
  Ytrain <- dat$Y[-test_set, ]
  Xtrain <- dat$X[-test_set, ]
  Ytest <- dat$Y[test_set, ]
  Xtest <- dat$X[test_set, ]
  
  return(list(Xtrain=Xtrain, Ytrain=Ytrain, Xtest=Xtest, Ytest=Ytest))
}