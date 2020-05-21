fit_mglmnet <- function(X, Y, alpha){
  time1 <- proc.time()
  
  cvfit <- cv.glmnet(x=X, y=Y, family="mgaussian", alpha=alpha)
  fit <- glmnet(x=X, y=Y, family="mgaussian", alpha=alpha, lambda=cvfit$lambda.min)
  
  time2 <- proc.time()
  
  Bhat <- matrix(as.numeric(NA), nrow=ncol(X), ncol=ncol(Y))
  for(i in 1:length(fit$beta)){
    Bhat[, i] <- as.vector(fit$beta[[i]])
  }
  
  elapsed_time <- time2["elapsed"] - time1["elapsed"]
  
  return(list(fit=fit, B_est=Bhat, intercept_est=drop(fit$a0), elapsed_time=elapsed_time))
}