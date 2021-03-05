###Function to compute initial estimates of the coefficients from group-lasso
compute_coefficients_glasso <- function(X, Y, standardize, nthreads, version=c("Rcpp", "R")){
  
  version <- match.arg(version)
  
  n <- nrow(X)
  p <- ncol(X)
  r <- ncol(Y)
  Y_has_missing <- any(is.na(Y))
  
  if(Y_has_missing){
    ###Extract per-individual Y missingness patterns
    Y_miss_patterns <- mr.mash.alpha:::extract_missing_Y_pattern(Y)
    
    ###Initialize missing Ys
    muy <- colMeans(Y, na.rm=TRUE)
    for(l in 1:r){
      Y[is.na(Y[, l]), l] <- muy[l]
    }
    
    ###Compute V and its inverse
    V <- mr.mash.alpha:::compute_V_init(X, Y, matrix(0, p, r), method="flash")
    Vinv <- chol2inv(chol(V))
    
    ###Compute expected Y (assuming B=0)
    mu <- matrix(rep(muy, each=n), n, r)
    
    ###Impute missing Ys 
    Y <- mr.mash.alpha:::impute_missing_Y(Y=Y, mu=mu, Vinv=Vinv, miss=Y_miss_patterns$miss, non_miss=Y_miss_patterns$non_miss, 
                                          version=version)$Y
  }
  
  ##Fit group-lasso
  if(nthreads>1){
    doMC::registerDoMC(nthreads)
    paral <- TRUE
  } else {
    paral <- FALSE
  }
  
  cvfit_glmnet <- glmnet::cv.glmnet(x=X, y=Y, family="mgaussian", alpha=1, standardize=standardize, parallel=paral)
  coeff_glmnet <- coef(cvfit_glmnet, s="lambda.min")
  
  ##Build matrix of initial estimates for mr.mash
  B <- matrix(as.numeric(NA), nrow=p, ncol=r)
  
  for(i in 1:length(coeff_glmnet)){
    B[, i] <- as.vector(coeff_glmnet[[i]])[-1]
  }
  
  return(list(B=B, Y=Y))
}
