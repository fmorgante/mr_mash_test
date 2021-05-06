simulate_data <- function(n, p, p_causal, r, r_causal,
                          pve, B_cor, B_scale, w, X_cor, 
                          X_scale, V_cor, prop_testset){

  if(!is.numeric(p_causal)){
    p_causal <- sample(x=eval(parse(text=p_causal)), size=1)
  }

  dat <- mr.mash.alpha::simulate_mr_mash_data(n=n, p=p, p_causal=p_causal, r=r,
                                               r_causal=r_causal, intercepts=rep(1, r),
                                               pve=pve, B_cor=B_cor, B_scale=B_scale,
                                               w=w, X_cor=X_cor, X_scale=X_scale,
                                               V_cor=V_cor)
  
  colnames(dat$Y) <- paste0("Y", seq(1, r))
  rownames(dat$Y) <- paste0("N", seq(1, n))
  colnames(dat$X) <- paste0("X", seq(1, p))
  rownames(dat$X) <- paste0("N", seq(1, n))
  
  if(prop_testset!=0){
    test_set <- sort(sample(x=c(1:n), size=round(n*prop_testset), replace=FALSE))
    Ytrain <- dat$Y[-test_set, ]
    Xtrain <- dat$X[-test_set, ]
    Ytest <- dat$Y[test_set, ]
    Xtest <- dat$X[test_set, ]
  
   return(list(Xtrain=Xtrain, Ytrain=Ytrain, Xtest=Xtest, Ytest=Ytest, B_true=dat$B))
  } else {
  return(list(X=dat$X, Y=dat$Y, B_true=dat$B))
  }
}
