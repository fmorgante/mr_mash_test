fit_mr_mash <- function(X, Y, update_w0, update_w0_method, standardize, update_V, ca_update_order, init_method, B_true,
                        w0_threshold, convergence_criterion, tol, daarem){
  
  r <- ncol(Y)
  p <- ncol(X)
  
  scaling_grid <- seq(0.1, 2.1, 0.2)


  ##Compute prior grid and weights for mr.ash (even if they will not be used)
  scaling_grid_mr_ash <- c(1e-10, scaling_grid)
  w0_mr_ash <- rep(1/length(scaling_grid_mr_ash), times=length(scaling_grid_mr_ash))
    
  if(init_method=="shared"){
    ##Prepare stacked genotypes and phenotypes
    y_mr_ash <- c()
    X_mr_ash <- X
    for(i in 1:r){
      y_mr_ash <- c(y_mr_ash, Y[, i])
      if(i != r)
        X_mr_ash <- rbind(X_mr_ash, X)
    }
      
    ##Fit mr.ash
    fit_mr_ash <- mr.ash.alpha::mr.ash(X_mr_ash, y_mr_ash, sa2=scaling_grid_mr_ash, pi=w0_mr_ash, 
                                       update.sigma=update_V, update.pi=update_w0, standardize=standardize, 
                                       verbose=FALSE, max.iter=5000)
      
    ##Build matrix of initial estimates for mr.mash
    mu1_init <- matrix(drop(fit_mr_ash$beta), nrow=p, ncol=r)
      
  } else if(init_method=="independent"){
    
    mu1_init <- matrix(as.numeric(NA), nrow=p, ncol=r)
    for(i in 1:r){
      y_mr_ash <- Y[, i]
        
      ##Fit mr.ash
      fit_mr_ash <- mr.ash.alpha::mr.ash(X, y_mr_ash, sa2=scaling_grid_mr_ash, pi=w0_mr_ash, 
                                         update.sigma=update_V, update.pi=update_w0, standardize=standardize, 
                                         verbose=FALSE, max.iter=5000)
        
      ##Build matrix of initial estimates for mr.mash
      mu1_init[, i] <- drop(fit_mr_ash$beta)
    }
  } else if(init_method=="2pass"){
      ##Prepare stacked genotypes and phenotypes
      y_mr_ash1 <- c()
      X_mr_ash1 <- X
      for(i in 1:r){
        y_mr_ash1 <- c(y_mr_ash1, Y[, i])
        if(i != r)
          X_mr_ash1 <- rbind(X_mr_ash1, X)
      }
    
      ##Fit mr.ash assuming shared effects
      fit_mr_ash1 <- mr.ash.alpha::mr.ash(X_mr_ash1, y_mr_ash1, sa2=scaling_grid_mr_ash, pi=w0_mr_ash, 
                                         update.sigma=update_V, update.pi=update_w0, standardize=standardize, 
                                         verbose=FALSE, max.iter=5000)
    
      ##Save initial estimates for mr.ash assuming independent effects
      beta_init1 <- drop(fit_mr_ash1$beta)
    
      ##Fit mr.ash assuming independent effects
      mu1_init <- matrix(as.numeric(NA), nrow=p, ncol=r)
      for(i in 1:r){
        y_mr_ash <- Y[, i]
      
        ##Fit mr.ash
        fit_mr_ash <- mr.ash.alpha::mr.ash(X, y_mr_ash, sa2=scaling_grid_mr_ash, pi=w0_mr_ash, 
                                           update.sigma=update_V, update.pi=update_w0, standardize=standardize, 
                                           verbose=FALSE, max.iter=5000, beta.init = beta_init1)
      
        ##Build matrix of initial estimates for mr.mash
        mu1_init[, i] <- drop(fit_mr_ash$beta)
      }
  } else if(init_method=="truth"){
    mu1_init <- B_true
  } else if(init_method=="default"){
    mu1_init <- matrix(0, nrow=p, ncol=r)
  } else if(init_method=="mlasso"){
    ##Fit multivariate LASSO
    cvfit_glmnet <- glmnet::cv.glmnet(x=X, y=Y, family="mgaussian", alpha=1)
    coeff_glmnet <- coef(cvfit_glmnet, s="lambda.min")
    
    ##Build matrix of initial estimates for mr.mash
    mu1_init <- matrix(as.numeric(NA), nrow=p, ncol=r)
    for(i in 1:length(coeff_glmnet)){
      mu1_init[, i] <- as.vector(coeff_glmnet[[i]])[-1]
    }
  }
  
  # if(init_method %in% c("shared", "independent", "2pass")){
  #   w0_up <- cbind(drop(fit_mr_ash$pi), scaling_grid_mr_ash)
  #   w0_up_sel <- w0_up[which(w0_up[, 1]>=select_w0_threshold), 2]
  #   scaling_grid <- w0_up_sel[-1] ##remove 0 because compute_cov_canonical() already has the null matrix
  # }
  
  ###Fit mr.mash
  S0 <- mr.mash.alpha:::compute_cov_canonical(r, singletons=TRUE, hetgrid=c(0, 0.25, 0.5, 0.75, 0.99), grid=scaling_grid, zeromat=TRUE)
  
  time1 <- proc.time()
  
  if(daarem){
    fit <- mr.mash.alpha::mr.mash.daar(X=X, Y=Y, S0=S0, update_w0=update_w0, update_w0_method=update_w0_method, 
                                       compute_ELBO=TRUE, standardize=standardize, verbose=FALSE, update_V=update_V, 
                                       version="Rcpp", ca_update_order=ca_update_order, mu1_init=mu1_init,
                                       mon_tol=1e-3, kappa=25, alpha=1.2, tol=tol)
  } else{
    fit <- mr.mash.alpha::mr.mash(X=X, Y=Y, S0=S0, update_w0=update_w0, update_w0_method=update_w0_method,
                                  compute_ELBO=TRUE, standardize=standardize, verbose=FALSE, update_V=update_V,
                                  version="Rcpp", ca_update_order=ca_update_order, mu1_init=mu1_init,
                                  w0_threshold=w0_threshold, convergence_criterion=convergence_criterion,
                                  tol=tol)
  }
  
  time2 <- proc.time()
  
  elapsed_time <- time2["elapsed"] - time1["elapsed"]
  
  return(list(fit=fit, B_est=fit$mu1, intercept_est=fit$intercept, elapsed_time=elapsed_time))
}
