fit_mr_mash_all_genes_prior <- function(X, Y, update_w0, update_w0_method, standardize, update_V, update_V_method, ca_update_order,
                                        w0_threshold, convergence_criterion, tol, data_driven_mats, sumstats, nthreads){
  
  ###Get number of responses and number of variables
  r <- ncol(Y)
  p <- ncol(X)
  
  time1 <- proc.time()
  
  ##Fit group-lasso
  cvfit_glmnet <- glmnet::cv.glmnet(x=X, y=Y, family="mgaussian", alpha=1, standardize=standardize)
  coeff_glmnet <- coef(cvfit_glmnet, s="lambda.min")
    
  ##Build matrix of initial estimates for mr.mash
  mu1_init <- matrix(as.numeric(NA), nrow=p, ncol=r)
  intercept <- rep(as.numeric(NA), r)
  for(i in 1:length(coeff_glmnet)){
    mu1_init[, i] <- as.vector(coeff_glmnet[[i]])[-1]
    intercept[i] <- as.vector(coeff_glmnet[[i]])[1]
  }
  
  ##Compute grid of scaling factors
  grid <- mr.mash.alpha::autoselect.mixsd(sumstats, mult=sqrt(2))^2
  
  if(is.null(data_driven_mats)){
    ###Compute canonical matrices
    S0_raw <- mr.mash.alpha::compute_canonical_covs(r, singletons=TRUE, hetgrid=c(0, 0.25, 0.5, 0.75, 1))
  } else {
    S0_raw <- readRDS(data_driven_mats)
  }
  
  ###Compute prior covariance and weights
  S0 <- mr.mash.alpha::expand_covs(S0_raw, grid, zeromat=TRUE)
  
  prop_nonzero_glmnet <- sum(mu1_init[, 1]!=0)/p
  w0 <- c((1-prop_nonzero_glmnet), rep(prop_nonzero_glmnet/(length(S0)-1), (length(S0)-1)))

  ###Fit mr.mash
  fit <- mr.mash.alpha::mr.mash(X=X, Y=Y, S0=S0, w0=w0, update_w0=update_w0, update_w0_method=update_w0_method,
                                compute_ELBO=TRUE, standardize=standardize, verbose=FALSE, update_V=update_V,
                                update_V_method=update_V_method, ca_update_order=ca_update_order, mu1_init=mu1_init, 
                                w0_threshold=w0_threshold, convergence_criterion=convergence_criterion,
                                tol=tol, nthreads=nthreads)
  
  time2 <- proc.time()
  
  elapsed_time <- time2["elapsed"] - time1["elapsed"]
  
  return(list(fit=fit, B_est=fit$mu1, intercept_est=fit$intercept, elapsed_time=elapsed_time))
}
