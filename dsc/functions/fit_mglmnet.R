fit_mglmnet <- function(X, Y, alpha){
  time1 <- proc.time()
  
  cvfit <- glmnet::cv.glmnet(x=X, y=Y, family="mgaussian", alpha=alpha)
  fit <- glmnet::glmnet(x=X, y=Y, family="mgaussian", alpha=alpha)
  coeffic <- coef(fit, s=cvfit$lambda.min)
  
  time2 <- proc.time()
  
  Bhat <- matrix(as.numeric(NA), nrow=(ncol(X)+1), ncol=ncol(Y))
  for(i in 1:length(coeffic)){
    Bhat[, i] <- as.vector(coeffic[[i]])
  }
  
  elapsed_time <- time2["elapsed"] - time1["elapsed"]
  
  return(list(fit=fit, B_est=Bhat[-1, ], intercept_est=Bhat[1, ], elapsed_time=elapsed_time))
}