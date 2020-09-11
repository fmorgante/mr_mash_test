simulate_data_all_genes_prior <- function(n, p, p_causal, r, r_causal,
                                          pve, B_cor, B_scale, w, X_cor, 
                                          X_scale, V_cor, testset_index){
  
  dat <- mr.mash.alpha::simulate_mr_mash_data(n=n, p=p, p_causal=p_causal, r=r,
                                               r_causal=r_causal, intercepts=rep(1, r),
                                               pve=pve, B_cor=B_cor, B_scale=B_scale,
                                               w=w, X_cor=X_cor, X_scale=X_scale,
                                               V_cor=V_cor)
  
  colnames(dat$Y) <- paste0("Y", seq(1, r))
  rownames(dat$Y) <- paste0("N", seq(1, n))
  colnames(dat$X) <- paste0("X", seq(1, p))
  rownames(dat$X) <- paste0("N", seq(1, n))
  
  test_set <- readRDS(testset_index)
  Ytrain <- dat$Y[-test_set, ]
  Xtrain <- dat$X[-test_set, ]
  Ytest <- dat$Y[test_set, ]
  Xtest <- dat$X[test_set, ]
  
  return(list(Xtrain=Xtrain, Ytrain=Ytrain, Xtest=Xtest, Ytest=Ytest, B_true=dat$B))
}