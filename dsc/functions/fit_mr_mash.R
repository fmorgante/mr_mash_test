fit_mr_mash <- function(X, Y, update_w0, update_w0_method, standardize, update_V, ca_update_order){
  
  r <- ncol(Y)
  S0 <- mr.mash.alpha:::compute_cov_canonical(r, singletons=TRUE, hetgrid=c(0, 0.25, 0.5, 0.75, 0.99), seq(0.1, 2, 0.2), zeromat=TRUE)
  
  time1 <- proc.time()
  fit <- mr.mash.alpha::mr.mash(X=X, Y=Y, S0=S0, update_w0=update_w0, update_w0_method=update_w0_method,
                 compute_ELBO=TRUE, standardize=standardize, verbose=FALSE, update_V=update_V,
                 version="Rcpp", ca_update_order=ca_update_order)
  time2 <- proc.time()
  elapsed_time <- time2["elapsed"] - time1["elapsed"]
  
  return(list(fit=fit, elapsed_time=elapsed_time))
}