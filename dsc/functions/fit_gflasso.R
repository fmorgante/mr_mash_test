fit_gflasso <- function(X, Y, nCores){
  X_s <- scale(X)
  Y_s <- scale(Y)
  
  time1 <- proc.time()
  
  cvfit <- gflasso::cv_gflasso(X = X_s, Y = Y_s, R = cor(Y)^2, nCores = nCores, seed=DSC_SEED)
  fit <- gflasso::gflasso(X = X_s, Y = Y_s, R = cor(Y)^2,
                         opts = list(lambda = cvfit$optimal$lambda, gamma = cvfit$optimal$gamma))
  
  time2 <- proc.time()
  
  Bhat <- t(t(fit$B)*attr(Y_s, 'scaled:scale'))/attr(X_s, 'scaled:scale')
  intercept <- attr(Y_s, 'scaled:center') - colSums(attr(X_s, 'scaled:center') * Bhat)
  
  elapsed_time <- time2["elapsed"] - time1["elapsed"]
  
  return(list(fit=fit, B_est=Bhat, intercept_est=intercept, elapsed_time=elapsed_time))
}