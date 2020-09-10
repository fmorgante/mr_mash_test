simulate_data <- function(n, p, p_causal, r, r_causal,
                          pve, B_cor, B_scale, w, X_cor, 
                          X_scale, V_cor, prop_testset){
  
  dat <- mr.mash.alpha::simulate_mr_mash_data(n=n, p=p, p_causal=p_causal, r=r,
                                               r_causal=r_causal, intercepts=rep(1, r),
                                               pve=pve, B_cor=B_cor, B_scale=B_scale,
                                               w=w, X_cor=X_cor, X_scale=X_scale,
                                               V_cor=V_cor)
  
  test_set <- sort(sample(x=c(1:n), size=round(n*prop_testset), replace=FALSE))
  Ytrain <- dat$Y[-test_set, ]
  Xtrain <- dat$X[-test_set, ]
  Ytest <- dat$Y[test_set, ]
  Xtest <- dat$X[test_set, ]
  
  return(list(Xtrain=Xtrain, Ytrain=Ytrain, Xtest=Xtest, Ytest=Ytest, B_true=dat$B))
}