###Load necessary libraries
library(mr.mash.alpha)
library(glmnet)
library(optparse)

options(stringsAsFactors = FALSE)

###Functions to compute accuracy and adjusted univariate sumstats
compute_accuracy <- function(Y, Yhat) {
  bias <- rep(as.numeric(NA), ncol(Y))
  r2 <- rep(as.numeric(NA), ncol(Y))
  mse <- rep(as.numeric(NA), ncol(Y))
  
  for(i in 1:ncol(Y)){
    fit  <- lm(Y[, i] ~ Yhat[, i])
    bias[i] <- coef(fit)[2] 
    r2[i] <- summary(fit)$r.squared
    mse[i] <- mean((Y[, i] - Yhat[, i])^2)
  }
  
  return(list(bias=bias, r2=r2, mse=mse))
}
compute_univariate_sumstats_adj <- function(X, Y, B, standardize=FALSE, standardize.response=FALSE){
  r <- ncol(Y)
  p <- ncol(X)
  Bhat <- matrix(as.numeric(NA), nrow=p, ncol=r)
  Shat <- matrix(as.numeric(NA), nrow=p, ncol=r)
  
  X <- scale(X, center=TRUE, scale=standardize) 
  Y <- scale(Y, center=TRUE, scale=standardize.response)
  
  if(standardize)
    B <- B*attr(X,"scaled:scale")
  
  for(i in 1:r){
    for(j in 1:p){
      Rij <- Y[, i] - X[, -j]%*%B[-j, i]
      fit <- lm(Rij ~ X[, j]-1)
      Bhat[j, i] <- coef(fit)
      Shat[j, i] <- summary(fit)$coefficients[1, 2]
    }
  }
  
  return(list(Bhat=Bhat, Shat=Shat))
}

###Parse command line arguments
parser <- OptionParser()
parser <- add_option(parser, c("--nthreads"), type="integer")
parser <- add_option(parser, c("--standardize"), type="logical")
parser <- add_option(parser, c("--w0_threshold"), type="numeric")
parser <- add_option(parser, c("--singletons"), type="logical")
outparse <- parse_args(parser)

###Set seed
set.seed(123)

###Set simulation parameters
n <- 900
p <- 5000
p_causal <- 5
r <- 50
r_causal <- list(1:10, 11:50)
B_cor <- c(1, 1)
B_scale <- c(1, 1)
w <- c(0.5, 0.5)

###Set some methods parameters
nthreads <- outparse$nthreads
w0_threshold <- outparse$w0_threshold
standardize <- outparse$standardize

###Simulate V, B, X and Y
out <- simulate_mr_mash_data(n, p, p_causal, r, r_causal, intercepts = rep(1, r),
                             pve=0.15, B_cor=B_cor, B_scale=B_scale, w=w,
                             X_cor=0.5, X_scale=1, V_cor=0)
colnames(out$Y) <- paste0("Y", seq(1, r))
rownames(out$Y) <- paste0("N", seq(1, n))
colnames(out$X) <- paste0("X", seq(1, p))
rownames(out$X) <- paste0("N", seq(1, n))

###Split the data in training and test sets
test_set <- sort(sample(x=c(1:n), size=round(n*0.2), replace=FALSE))
Ytrain <- out$Y[-test_set, ]
Xtrain <- out$X[-test_set, ]
Ytest <- out$Y[test_set, ]
Xtest <- out$X[test_set, ]

###Fit grop-lasso to initialize mr.mash
tic <- Sys.time()
cvfit_glmnet <- cv.glmnet(x=Xtrain, y=Ytrain, family="mgaussian", alpha=1, standardize=standardize)
toc <- Sys.time()

cat("group-lasso executed in", difftime(toc, tic, units="mins"),
    "minutes\n")

coeff_glmnet <- coef(cvfit_glmnet, s="lambda.min")
Bhat_glmnet <- matrix(as.numeric(NA), nrow=p, ncol=r)
for(i in 1:length(coeff_glmnet)){
  Bhat_glmnet[, i] <- as.vector(coeff_glmnet[[i]])[-1]
}
Yhat_glmnet <- drop(predict(cvfit_glmnet, newx=Xtest, s="lambda.min"))
prop_nonzero_glmnet <- sum(Bhat_glmnet[, 1]!=0)/p

###Compute grid of variances
tic <- Sys.time()
univ_sumstats <- compute_univariate_sumstats(Xtrain, Ytrain, standardize=standardize, standardize.response=FALSE)
grid <- autoselect.mixsd(univ_sumstats, mult=sqrt(2))^2
toc <- Sys.time()

cat("grid of variances computed in", difftime(toc, tic, units="mins"),
    "minutes\n")

###Compute prior with both canonical and data-driven covariance matrices
tic <- Sys.time()
univ_sumstats_adj <- compute_univariate_sumstats_adj(Xtrain, Ytrain, Bhat_glmnet, standardize=standardize, standardize.response=FALSE)
res_cor <- mashr::estimate_null_correlation_simple(mashr::mash_set_data(univ_sumstats_adj$Bhat, univ_sumstats_adj$Shat))
S0_datadriven <- compute_data_driven_covs(univ_sumstats_adj, subset_thresh=0.05, n_pcs=3, 
                                          non_singleton=FALSE, Gamma=res_cor)
S0_can <- compute_canonical_covs(ncol(Ytrain), singletons=outparse$singletons, hetgrid=c(0, 0.5, 1))
S0_full <- c(S0_can, S0_datadriven)
S0_full <- expand_covs(S0_full, grid, zeromat=TRUE)
toc <- Sys.time()

cat("prior computed in", difftime(toc, tic, units="mins"),
    "minutes\n")

###Fit mr.mash
w0_full <- c((1-prop_nonzero_glmnet), rep(prop_nonzero_glmnet/(length(S0_full)-1), (length(S0_full)-1)))
fit_mrmash_full <- mr.mash(Xtrain, Ytrain, S0_full, w0=w0_full, update_w0=TRUE, update_w0_method="EM", tol=1e-2,
                           convergence_criterion="ELBO", compute_ELBO=TRUE, standardize=standardize, 
                           verbose=TRUE, update_V=TRUE, update_V_method="full", e=1e-8,
                           mu1_init=Bhat_glmnet, w0_threshold=w0_threshold, nthreads=nthreads)
Yhat_mrmash_full <- predict(fit_mrmash_full, Xtest)

###Accuracy of prediction
cat("Mean scaled MSE for group-lasso", mean(compute_accuracy(Ytest, Yhat_glmnet)$mse/matrixStats::colVars(Ytest)), "\n")
cat("Mean scaled MSE for mr.mash", mean(compute_accuracy(Ytest, Yhat_mrmash_full)$mse/matrixStats::colVars(Ytest)), "\n")


