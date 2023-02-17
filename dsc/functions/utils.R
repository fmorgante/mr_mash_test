###Functions to compute adjusted univariate sumstats
compute_univariate_sumstats_adj <- function(X, Y, B, a, standardize=standardize, standardize.response=FALSE, mc.cores=1){
  r <- ncol(Y)
  
  X <- scale(X, center=TRUE, scale=standardize) 
  Y <- scale(t(t(Y)-a), center=TRUE, scale=standardize.response)
  
  if(standardize)
    B <- B*attr(X,"scaled:scale")
  
  linreg <- function(i, X, Y, B){
    p <- ncol(X)
    bhat <- rep(as.numeric(NA), p)
    shat <- rep(as.numeric(NA), p)
    
    for(j in 1:p){
      Rij <- Y[, i] - X[, -j]%*%B[-j, i]
      fit <- lm(Rij ~ X[, j]-1)
      bhat[j] <- coef(fit)
      shat[j] <- summary(fit)$coefficients[1, 2]
    }
    
    return(list(bhat=bhat, shat=shat))
  }
  
  ###mclapply is a little faster but uses more memory
  #out <- parallel::mclapply(1:r, linreg, X, Y, B, mc.cores=mc.cores)
  
  cl <- parallel::makeCluster(mc.cores)
  out <- parallel::parLapply(cl, 1:r, linreg, X, Y, B)
  parallel::stopCluster(cl)
  
  return(list(Bhat=sapply(out,"[[","bhat"), Shat=sapply(out,"[[","shat")))
}


###Functions to compute MAF and missing genotype rate
compute_maf <- function(geno){
  f <- mean(geno,na.rm = TRUE)/2
  return(min(f, 1-f))
}

compute_missing <- function(geno){
  miss <- sum(is.na(geno))/length(geno)
  return(miss)
}

mean_impute <- function(geno){
  f <- apply(geno, 2, function(x) mean(x,na.rm = TRUE))
  for (i in 1:length(f)) geno[,i][which(is.na(geno[,i]))] <- f[i]
  return(geno)
}

is_zero_variance <- function(x) {
  if (length(unique(x[!is.na(x)]))==1) return(T)
  else return(F)
}


### Filter X matrix
filter_X <- function(X, missing_rate_thresh, maf_thresh, var_thresh) {
  rm_col <- which(apply(X, 2, compute_missing) > missing_rate_thresh)
  if (length(rm_col)) X <- X[, -rm_col]
  rm_col <- which(apply(X, 2, compute_maf) < maf_thresh)
  if (length(rm_col)) X <- X[, -rm_col]
  X <- mean_impute(X)
  #rm_col <- which(apply(X, 2, is_zero_variance))
  rm_col <- which(matrixStats::colVars(X) < var_thresh)
  if (length(rm_col)) X <- X[, -rm_col]
  return(X)
}

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

###Function to assign missing values to Y randomly
assign_NAs_random <- function(Y, prop){
  tot_num <- length(Y)
  miss_num <- floor(tot_num*prop)
  miss_idx <- sort(sample(x=1:tot_num, size=miss_num))
  
  Yna <- t(Y)
  Yna[miss_idx] <- NA
  
  return(t(Yna))
}

###Function to impute missing values with the column mean
col_mean_impute <- function(mat){
  
  means <- colMeans(mat, na.rm=TRUE)
  r <- ncol(mat)
  
  for(i in 1:r){
    whichNa <- which(is.na(mat[,i]))
    mat[whichNa, i] <- means[i]
  }
  
  return(mat)
}

