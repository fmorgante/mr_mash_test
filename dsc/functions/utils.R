###Functions to compute adjusted univariate sumstats
compute_univariate_sumstats_adj <- function(X, Y, B, a, standardize=FALSE, standardize.response=FALSE, mc.cores=1){
  r <- ncol(Y)
  
  X <- scale(X, center=TRUE, scale=standardize) 
  Y <- scale(t(t(Y)-a), center=TRUE, scale=standardize.response)
  
  if(standardize)
    B <- B*attr(X,"scaled:scale")
  
  linreg <- function(i, X, Y, B){
    p <- ncol(X)
    bhat <- rep(as.numeric(NA), p)
    shat <- rep(as.numeric(NA), p)
    
    for(j in 1:p){
      Rij <- Y[, i] - X[, -j]%*%B[-j, i]
      fit <- lm(Rij ~ X[, j]-1)
      bhat[j] <- coef(fit)
      shat[j] <- summary(fit)$coefficients[1, 2]
    }
    
    return(list(bhat=bhat, shat=shat))
  }
  
  ###mclapply is a little faster but uses more memory
  #out <- parallel::mclapply(1:r, linreg, X, Y, B, mc.cores=mc.cores)
  
  cl <- parallel::makeCluster(mc.cores)
  out <- parallel::parLapply(cl, 1:r, linreg, X, Y, B)
  parallel::stopCluster(cl)
  
  return(list(Bhat=sapply(out,"[[","bhat"), Shat=sapply(out,"[[","shat")))
}
