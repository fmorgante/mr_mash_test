fit_utmost <- function(X, Y, nfolds){
  Y <- scale(Y, scale=FALSE)
  muy <- attr(Y, 'scaled:center')
  mux <- colMeans(X)

  time1 <- proc.time()

  B <- detect_eqtls_glasso(X, Y, nfolds)
  intercept <- drop(muy - mux %*% B)

  time2 <- proc.time()

  elapsed_time <- time2["elapsed"] - time1["elapsed"]

  return(list(B_est=B, intercept_est=intercept, elapsed_time=elapsed_time))
}
