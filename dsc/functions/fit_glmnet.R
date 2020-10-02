fit_glmnet <- function(X, Y, alpha, standardize, nthreads){
  
  r <- ncol(Y)
  
  linreg <- function(i, X, Y, alpha, standardize, nthreads){
    if(nthreads>1){
      doMC::registerDoMC(nthreads)
      paral <- TRUE
    } else {
      paral <- FALSE
    }
    
    cvfit <- glmnet::cv.glmnet(x=X, y=Y[, i], family="gaussian", alpha=alpha, standardize=standardize, parallel=paral)
    coeffic <- as.vector(coef(cvfit, s="lambda.min"))
    lambda_seq <- cvfit$lambda
    
    return(list(bhat=coeffic, lambda_seq=lambda_seq))
  }
  
  time1 <- proc.time()
  
  out <- lapply(1:r, linreg, X, Y, alpha, standardize, nthreads)
  
  time2 <- proc.time()

  Bhat <- sapply(out,"[[","bhat")
  lambda_seq <- unlist(lapply(out,"[[","lambda_seq"))
  lambda_maxmin <- c(max(lambda_seq), min(lambda_seq))
  elapsed_time <- time2["elapsed"] - time1["elapsed"]

  return(list(B_est=Bhat[-1, ], intercept_est=Bhat[1, ], lambda_maxmin=lambda_maxmin, elapsed_time=elapsed_time))
}
