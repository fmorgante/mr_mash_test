fit_mglmnet <- function(X, Y, alpha, standardize, nthreads){
  
  if(nthreads>1){
    doMC::registerDoMC(nthreads)
    paral <- TRUE
  } else {
    paral <- FALSE
  }
  
  time1 <- proc.time()
  
  cvfit <- glmnet::cv.glmnet(x=X, y=Y, family="mgaussian", alpha=alpha, standardize=standardize, parallel=paral)
  coeffic <- coef(cvfit, s="lambda.min")
  
  time2 <- proc.time()
  
  Bhat <- matrix(as.numeric(NA), nrow=(ncol(X)+1), ncol=ncol(Y))
  for(i in 1:length(coeffic)){
    Bhat[, i] <- as.vector(coeffic[[i]])
  }
  
  elapsed_time <- time2["elapsed"] - time1["elapsed"]
  
  return(list(fit=cvfit, B_est=Bhat[-1, ], intercept_est=Bhat[1, ], elapsed_time=elapsed_time))
}
