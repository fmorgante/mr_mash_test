set.seed(123)

###Load library
library(mr.mash.alpha)
#library(glmnet)

###Simulate data from given X
simulate_mr_mash_data_from_given_X <- function(X, p_causal, r, r_causal, intercepts,
                                               pve, B_cor, B_scale, w, V_cor){
  ##Check that the inputs are correct
  if(!is.matrix(X))
    stop("X must be a matrix.")
  if(any(is.na(X)))
    stop("X must not contain missing values.")
  if(length(intercepts)!=r)
    stop("intercepts must be of length equal to r.")
  if(any(sapply(r_causal, length)>r))
    stop("r_causal cannot be greater than r.")
  if(!(length(B_cor)==length(B_scale) & length(w)==length(B_cor) & length(B_scale)==length(w)))
    stop("B_cor, B_scale, and w must have the same length.")
  if(abs(sum(w) - 1) > 1e-10)
    stop("Elements of w must sum to 1.")
  if(length(pve)!=1 & length(pve)!=r)
    stop("pve must be of length equal to 1 or r.")
  
  ##Get number of mixture components, samples, and variables
  n <- nrow(X)
  p <- ncol(X)
  K <- length(w)
  
  ##Simulate true effects from N_r(0, Sigma) or \sum_K w_k N_r(0, Sigma_k) where Sigma and Sigma_k are given 
  ##covariance matrices across traits and w_k is the mixture proportion associated to Sigma_k
  Sigma <- vector("list", K)
  for(i in 1:K){
    r_mix_length <- length(r_causal[[i]])
    Sigma_offdiag <- B_scale[i]*B_cor[i]
    Sigma[[i]] <- matrix(Sigma_offdiag, nrow=r_mix_length, ncol=r_mix_length)
    diag(Sigma[[i]]) <- B_scale[i]
  }
  #Sample effects from a mixture of MVN distributions or a single MVN distribution
  B_causal <- matrix(0, nrow=p_causal, ncol=r)
  if(K>1){
    mixcomps <- sample(x=1:K, size=p_causal, prob=w, replace=TRUE)
    for(j in 1:p_causal){
      comp_to_use <- mixcomps[j]
      r_causal_mix <- r_causal[[comp_to_use]]
      B_causal[j, r_causal_mix] <- mvtnorm::rmvnorm(n=1, mean=rep(0, length(r_causal_mix)), sigma=Sigma[[comp_to_use]])
    }
  } else {
    r_causal_length <- length(r_causal[[1]])
    r_causal_index <- r_causal[[1]]
    B_causal[, r_causal_index] <- mvtnorm::rmvnorm(n=p_causal, mean=rep(0, r_causal_length), sigma=Sigma[[1]])
  }
  B <- matrix(0, ncol=r, nrow=p)
  causal_variables <- sample(x=(1:p), size=p_causal)
  B[causal_variables, ] <- B_causal
  
  ##Center X
  Xc <- mr.mash.alpha:::scale_fast2(X, scale=FALSE)$M
  
  ##Compute G and its variance
  G <- Xc%*%B
  Var_G <- matrixStats::colVars(G)
  
  ##Compute residual covariance
  Var_E <- ((1/pve)-1)*Var_G
  Var_E[which(Var_E<=.Machine$double.eps)] <- 1
  D <- diag(x=sqrt(Var_E))
  V_cor_mat <- matrix(V_cor, nrow=r, ncol=r)
  diag(V_cor_mat) <- 1
  V <- D %*% V_cor_mat %*% D
  
  ##Simulate Y from MN(XB, I_n, V) where I_n is an nxn identity matrix and V is the residual covariance  
  Y <- MBSP::matrix.normal(G + matrix(intercepts, n, r, byrow=TRUE), diag(n), V)
  
  ##Compile output
  causal_responses <- r_causal
  names(causal_responses) <- paste0("Component", 1:K)
  names(Sigma) <- paste0("Component", 1:K)
  out <- list(X=Xc, X_orig=X, Y=Y, B=B, V=V, Sigma=Sigma, intercepts=intercepts, causal_responses=causal_responses)
  if(K>1){
    if(p_causal>1){
      causal_variables_mixcomps <- cbind(causal_variables, mixcomps)
      causal_variables_mixcomps <- causal_variables_mixcomps[order(causal_variables_mixcomps[, 1]), ]
      out$causal_variables <- causal_variables_mixcomps[, 1]
      out$causal_vars_to_mixture_comps <- causal_variables_mixcomps[, 2]
    } else {
      out$causal_variables <- causal_variables
      out$causal_vars_to_mixture_comps <- mixcomps
    }
  } else {
    out$causal_variables <- sort(causal_variables)
    out$causal_vars_to_mixture_comps <- rep(1, p_causal)
  }
  
  return(out)
}

fit_glmnet <- function(X, Y, alpha, standardize, nthreads){
  
  r <- ncol(Y)
  
  linreg <- function(i, X, Y, alpha, standardize, nthreads){
    if(nthreads>1){
      doMC::registerDoMC(nthreads)
      paral <- TRUE
    } else {
      paral <- FALSE
    }
    
    cvfit <- glmnet::cv.glmnet(x=X, y=Y[, i], family="gaussian", alpha=alpha, standardize=standardize, parallel=paral)
    coeffic <- as.vector(coef(cvfit, s="lambda.min"))
    lambda_seq <- cvfit$lambda
    
    return(list(bhat=coeffic, lambda_seq=lambda_seq))
  }
  
  time1 <- proc.time()
  
  out <- lapply(1:r, linreg, X, Y, alpha, standardize, nthreads)
  
  time2 <- proc.time()

  Bhat <- sapply(out,"[[","bhat")
  lambda_seq <- unlist(lapply(out,"[[","lambda_seq"))
  lambda_maxmin <- c(max(lambda_seq), min(lambda_seq))
  elapsed_time <- time2["elapsed"] - time1["elapsed"]

  return(list(B_est=Bhat[-1, ], intercept_est=Bhat[1, ], lambda_maxmin=lambda_maxmin, elapsed_time=elapsed_time))
}

###Load wtccc data
load("../data/t1d.RData")

###Simulate data
dat <- simulate_mr_mash_data_from_given_X(X=X, p_causal=1000, r=10, r_causal=list(1:10), 
                                          intercepts=rep(1, 10), pve=0.5, B_cor=0.8, 
                                          B_scale=1, w=1, V_cor=0)
Y <- dat$Y
B <- dat$B
V <- dat$V

rm(list=c("chr", "labels", "major", "minor", "pos", "y", "dat"))
gc()

###Compute initial B
#B_glmnet <- fit_glmnet(X=X, Y=Y, alpha=0.5, standardize=TRUE, nthreads=4)$B_est
gc()

###Compute prior
univ_sumstats <- compute_univariate_sumstats(X=X, Y=Y, standardize=TRUE, 
                                             standardize.response=FALSE, mc.cores=1)
grid <- autoselect.mixsd(univ_sumstats, mult=sqrt(2))^2
S0_can <- compute_canonical_covs(ncol(Y), singletons=TRUE, hetgrid=c(0, 0.25, 0.5, 0.75, 1))
S0 <- expand_covs(S0_can, grid, zeromat=TRUE)

rm(univ_sumstats)

###Fit mr.mash
w0 <- c(0.9, rep(0.1/(length(S0)-1), (length(S0))-1))

fit_mrmash <- mr.mash(X=X, Y=Y, S0=S0, w0=w0, update_w0=TRUE, tol=1e-2,
                      convergence_criterion="ELBO", compute_ELBO=TRUE, standardize=TRUE, 
                      verbose=TRUE, update_V=TRUE, update_V_method="diagonal", e=1e-8,
                      w0_threshold=1e-8, nthreads=8)
                      #, mu1_init=B_glmnet)

###Save results
res <- list(X=X, Y=Y, B=B, V=V, fit_mrmash=fit_mrmash) #, B_glmnet=B_glmnet)

saveRDS(res, "Desktop/res_wtccc_mrmash_0.rds")

