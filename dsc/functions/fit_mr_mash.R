fit_mr_mash <- function(X, Y, update_w0, update_w0_method, standardize, update_V, ca_update_order, mr_ash_method){
  
  r <- ncol(Y)
  p <- ncol(X)
  
  scaling_grid <- seq(0.1, 2.1, 0.2)

  ###Fit mr.ash to get initial estimates of mu1 for mr.mash, if requested
  if(!is.null(mr_ash_method)){
    ##Compute prior grid and weights for mr.ash
    scaling_grid_mr_ash <- c(1e-10, scaling_grid)
    w0_mr_ash <- rep(1/length(scaling_grid_mr_ash), times=length(scaling_grid_mr_ash))
    
    if(mr_ash_method=="shared"){
      ##Prepare stacked genotypes and phenotypes
      y_mr_ash <- c()
      X_mr_ash <- X
      for(i in 1:r){
        y_mr_ash <- c(y_mr_ash, Y[, i])
        if(i != r)
          X_mr_ash <- rbind(X_mr_ash, X)
      }
      
      ##Fit mr.ash
      fit_mr_ash <- mr.ash.alpha::mr.ash(X_mr_ash, y_mr_ash, sigma2=var(y_mr_ash), sa2=scaling_grid_mr_ash/var(y_mr_ash), pi=w0_mr_ash, 
                                         update.sigma=update_V, update.pi=update_w0, standardize=standardize, verbose=FALSE, max.iter=5000)
      
      ##Build matrix of initial estimates for mr.mash
      mu1_init <- matrix(fit_mr_ash$beta, nrow=p, ncol=r)
      
    } else if(mr_ash_method=="independent"){
      
      mu1_init <- matrix(as.numeric(NA), nrow=p, ncol=r)
      for(i in 1:r){
        y_mr_ash <- Y[, i]
        
        ##Fit mr.ash
        fit_mr_ash <- mr.ash.alpha::mr.ash(X, y_mr_ash, sigma2=var(y_mr_ash), sa2=scaling_grid_mr_ash/var(y_mr_ash), pi=w0_mr_ash, 
                                           update.sigma=update_V, update.pi=update_w0, standardize=standardize, verbose=FALSE, max.iter=5000)
        
        ##Build matrix of initial estimates for mr.mash
        mu1_init[, i] <- fit_mr_ash$beta
      }
    }
  } else {
    mu1_init <- matrix(0, nrow=p, ncol=r)
  }
  
  ###Fit mr.mash
  S0 <- mr.mash.alpha:::compute_cov_canonical(r, singletons=TRUE, hetgrid=c(0, 0.25, 0.5, 0.75, 0.99), grid=scaling_grid, zeromat=TRUE)
  
  time1 <- proc.time()
  
  fit <- mr.mash.alpha::mr.mash(X=X, Y=Y, S0=S0, update_w0=update_w0, update_w0_method=update_w0_method,
                                  compute_ELBO=TRUE, standardize=standardize, verbose=FALSE, update_V=update_V,
                                  version="Rcpp", ca_update_order=ca_update_order, mu1_init=mu1_init)
  
  time2 <- proc.time()
  
  elapsed_time <- time2["elapsed"] - time1["elapsed"]
  
  return(list(fit=fit, B_est=fit$mu1, intercept_est=fit$intercept, elapsed_time=elapsed_time))
}
