fit_utmost <- function(X, Y, nfolds){
  time1 <- proc.time()

  B <- detect_eqtls_glasso(X, Y, nfolds)

  time2 <- proc.time()

  elapsed_time <- time2["elapsed"] - time1["elapsed"]

  return(list(B_est=B, elapsed_time=elapsed_time))
}
