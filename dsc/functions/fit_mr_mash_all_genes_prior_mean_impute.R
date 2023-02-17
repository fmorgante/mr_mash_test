fit_mr_mash_all_genes_prior <- function(X, Y, update_w0, update_w0_method, standardize, update_V, update_V_method, ca_update_order,
                                        w0_threshold, convergence_criterion, tol, canonical_mats, singletons, hetgrid, data_driven_mats, sumstats, 
                                        nthreads, mu1_init){
  
  ###Get number of responses and number of variables
  r <- ncol(Y)
  p <- ncol(X)
  
  time1 <- proc.time()
  
  ##Compute grid of scaling factors
  grid <- mr.mash.alpha::autoselect.mixsd(sumstats, mult=sqrt(2))^2
  
  ###Compute canonical matrices
  if(canonical_mats){
    S0_raw <- mr.mash.alpha::compute_canonical_covs(r, singletons=singletons, hetgrid=hetgrid)
  }
  
  ###Load and extract data-driven matrices, if requested
  if(!is.null(data_driven_mats)){
    S0_data <- readRDS(data_driven_mats)
    S0_data[c("identity", colnames(Y), "equal_effects", "simple_het_1", "simple_het_2", "simple_het_3")] <- NULL
    if(canonical_mats){
      S0_raw <- c(S0_raw, S0_data)
    } else {
      S0_raw <- S0_data
    }
  }
  
  ###Compute prior covariance and weights
  S0 <- mr.mash.alpha::expand_covs(S0_raw, grid, zeromat=TRUE)
  
  prop_nonzero_glmnet <- sum(rowSums(abs(mu1_init))>0)/p
  w0 <- c((1-prop_nonzero_glmnet), rep(prop_nonzero_glmnet/(length(S0)-1), (length(S0)-1)))

  ###Mean impute missing Y
  Y <- col_mean_impute(Y)

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
