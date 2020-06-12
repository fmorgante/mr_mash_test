library(mr.mash.alpha)
library(glmnet)

setwd("/project2/mstephens/fmorgante/mr_mash_test")

###Set options
options(stringsAsFactors = FALSE)

###Function to compute accuracy
compute_accuracy <- function(Y, Yhat) {
  if(!is.null(dim(Y))){
    bias <- rep(as.numeric(NA), ncol(Y))
    r2 <- rep(as.numeric(NA), ncol(Y))
    mse <- rep(as.numeric(NA), ncol(Y))
    
    for(i in 1:ncol(Y)){
      fit  <- lm(Y[, i] ~ Yhat[, i])
      bias[i] <- coef(fit)[2] 
      r2[i] <- summary(fit)$r.squared
      mse[i] <- mean((Y[, i] - Yhat[, i])^2)
    }
  } else {
    fit  <- lm(Y ~ Yhat)
    bias <- coef(fit)[2] 
    r2 <- summary(fit)$r.squared
    mse <- mean((Y - Yhat)^2)
  }
  
  return(list(bias=bias, r2=r2, mse=mse))
}

###Create quantities to stores the results
lasso_mse <- rep(as.numeric(NA), 10)
group_lasso_mse <- matrix(as.numeric(NA), 10, 5)
mr_mash_FE <- matrix(as.numeric(NA), 10, 5)
mr_mash_FE_diagV <- matrix(as.numeric(NA), 10, 5)
mr_ash_mse <- rep(as.numeric(NA), 10)
mr_mash_1d_mse <- rep(as.numeric(NA), 10)

z <- 0

###Loop over seeds
RNGversion("3.5.0")
for(seed in c(1, 65, 23, 5009, 2857, 33, 579, 319, 498, 1040)){
  ###Set seed
  set.seed(seed)
  
  ###Indices needed to place the results in the correct position  
  z <- z+1

  ###Set parameters
  n <- 1000
  p <- 1000
  p_causal <- 500
  r <- 5
  r_causal <- r
  pve <- 0.5
  B_cor <- 1
  X_cor <- 0
  V_cor <- 0
  
  ###Simulate V, B, X and Y
  out <- mr.mash.alpha:::simulate_mr_mash_data(n, p, p_causal, r, r_causal=r, intercepts = rep(1, r),
                                               pve=pve, B_cor=B_cor, B_scale=0.9, X_cor=X_cor, X_scale=0.8, V_cor=V_cor)
  
  colnames(out$Y) <- paste0("Y", seq(1, r))
  rownames(out$Y) <- paste0("N", seq(1, n))
  colnames(out$X) <- paste0("X", seq(1, p))
  rownames(out$X) <- paste0("N", seq(1, n))
  
  ###Split the data in training and test sets
  test_set <- sort(sample(x=c(1:n), size=round(n*0.5), replace=FALSE))
  Ytrain <- out$Y[-test_set, ]
  Xtrain <- out$X[-test_set, ]
  Ytest <- out$Y[test_set, ]
  Xtest <- out$X[test_set, ]
  
  ###Fit group-lasso
  cvfit_glmnet <- cv.glmnet(x=Xtrain, y=Ytrain, family="mgaussian", alpha=1, standardize=FALSE)
  coeff_glmnet <- coef(cvfit_glmnet, s="lambda.min")
  Bhat_glmnet <- matrix(as.numeric(NA), nrow=p, ncol=r)
  ahat_glmnet <- vector("numeric", r)
  for(i in 1:length(coeff_glmnet)){
    Bhat_glmnet[, i] <- as.vector(coeff_glmnet[[i]])[-1]
    ahat_glmnet[i] <- as.vector(coeff_glmnet[[i]])[1]
  }

  ###Compute row means of Y
  Ymean <- cbind(rowSums(Ytrain)/ncol(Ytrain))
  
  ###Fit lasso on row the row means
  cvfit_glmnet_oned <- cv.glmnet(x=Xtrain, y=Ymean, family="gaussian", alpha=1, standardize=FALSE)
  coeff_glmnet_oned <- coef(cvfit_glmnet_oned, s="lambda.min")
  bhat_glmnet_oned <- as.vector(coeff_glmnet_oned)[-1]
  ahat_glmnet_oned <- as.vector(coeff_glmnet_oned)[1]

  ###Fit mr.mash.FE
  univ_sumstats <- mr.mash.alpha:::get_univariate_sumstats(Xtrain, Ytrain, standardize=FALSE, standardize.response=FALSE)
  grid <- mr.mash.alpha:::autoselect.mixsd(c(univ_sumstats$Bhat), c(univ_sumstats$Shat), mult=sqrt(2))^2
  S0 <- mr.mash.alpha:::compute_cov_canonical(ncol(Ytrain), singletons=FALSE, hetgrid=c(1), grid, zeromat=TRUE)
  fit_mrmashFE <- mr.mash(Xtrain, Ytrain, S0, tol=1e-2, convergence_criterion="ELBO", update_w0=TRUE,
                          update_w0_method="EM", compute_ELBO=TRUE, standardize=FALSE, verbose=FALSE, update_V=TRUE,
                          mu1_init=matrix(bhat_glmnet_oned, nrow=p, ncol=r), w0_threshold=1e-8)

  ###Fit mr.mash.FE with diagonal V
  univ_sumstats <- mr.mash.alpha:::get_univariate_sumstats(Xtrain, Ytrain, standardize=FALSE, standardize.response=FALSE)
  grid <- mr.mash.alpha:::autoselect.mixsd(c(univ_sumstats$Bhat), c(univ_sumstats$Shat), mult=sqrt(2))^2
  S0 <- mr.mash.alpha:::compute_cov_canonical(ncol(Ytrain), singletons=FALSE, hetgrid=c(1), grid, zeromat=TRUE)
  fit_mrmashFE_diagV <- mr.mash(Xtrain, Ytrain, S0, tol=1e-2, convergence_criterion="ELBO", update_w0=TRUE,
                                update_w0_method="EM", compute_ELBO=TRUE, standardize=FALSE, verbose=FALSE, update_V=TRUE,
                                update_V_method="diagonal", mu1_init=matrix(bhat_glmnet_oned, nrow=p, ncol=r), w0_threshold=1e-8)
  
  ###Fit mr.ash
  fit_mrash <- mr.ash.alpha::mr.ash(Xtrain, Ymean, beta.init=bhat_glmnet_oned) 
  
  ###Fit mr.mash univariate
  s0 <- vector("list", length(grid)+1)
  for(i in 1:(length(grid)+1)){
    s0[[i]] <- matrix(c(0, grid)[i], ncol=1, nrow=1)
  }
  fit_mrmash_1d <- mr.mash(Xtrain, Ymean, s0, tol=1e-2, convergence_criterion="ELBO", update_w0=TRUE,
                          update_w0_method="EM", compute_ELBO=TRUE, standardize=FALSE, verbose=FALSE, update_V=TRUE,
                          mu1_init=matrix(bhat_glmnet_oned, nrow=p, ncol=1), w0_threshold=0)
  
  ###Compute MSE
  group_lasso_mse[z, ] <- compute_accuracy(out$B, Bhat_glmnet)$mse
  lasso_mse[z] <- compute_accuracy(out$B[,1], bhat_glmnet_oned)$mse
  mr_mash_FE[z, ] <- compute_accuracy(out$B, fit_mrmashFE$mu1)$mse
  mr_mash_FE_diagV[z, ] <- compute_accuracy(out$B, fit_mrmashFE_diagV$mu1)$mse
  mr_ash_mse[z] <- compute_accuracy(out$B[,1], fit_mrash$beta)$mse
  mr_mash_1d_mse[z] <- compute_accuracy(out$B[,1], fit_mrmash_1d$mu1)$mse
}


group_lasso_mse

mr_mash_FE

mr_mash_FE_diagV

lasso_mse

mr_ash_mse

mr_mash_1d_mse
