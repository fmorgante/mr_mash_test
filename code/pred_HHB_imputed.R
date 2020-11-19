options(stringsAsFactors = FALSE)

library(mr.mash.alpha)
library(glmnet)

set.seed(1)

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

###Functions to compute accuracy
compute_accuracy <- function(Y, Yhat) {
  bias <- rep(as.numeric(NA), ncol(Y))
  r2 <- rep(as.numeric(NA), ncol(Y))
  rmse <- rep(as.numeric(NA), ncol(Y))
  
  for(i in 1:ncol(Y)){
    if(var(Yhat[, i])>0){
      fit  <- lm(Y[, i] ~ Yhat[, i])
      bias[i] <- coef(fit)[2] 
      r2[i] <- summary(fit)$r.squared
      rmse[i] <- sqrt(mean((Y[, i] - Yhat[, i])^2))
    } else {
      bias[i] <- NA
      r2[i] <- NA
      rmse[i] <- sqrt(mean((Y[, i] - Yhat[, i])^2))
    }
    
  }
  
  return(list(bias=bias, r2=r2, rmse=rmse, scaled_rmse=rmse/matrixStats::colSds(Y)))
}

###Load data 
dat <- readRDS("/project2/mstephens/fmorgante/mr_mash_test/data/cis_eqtl_analysis_ready/ENSG00000244734.3.Multi_Tissues.rds")
Y <- dat$y_res

###Impute missing values with flash
apply(Y, 2, function(x) all(is.na(x)))
fits <- flashier::flash(Y)
Y[is.na(Y)] <- fitted(fits)[is.na(Y)]

###Filter X
X <- filter_X(dat$X, 0.05, 0.05, 0.05)

###Split the data in training and test sets
test_set <- sort(sample(x=rownames(X), size=168, replace=FALSE))

Ytrain <- Y[!(rownames(Y) %in% test_set), ]
Xtrain <- X[!(rownames(X) %in% test_set), ]
Ytest <- Y[test_set, ]
Xtest <- X[test_set, ]

###Set some parameters
standardize <- TRUE
nthreads <- 8
w0_threshold <- 1e-8

###Compute grid of variances
univ_sumstats <- compute_univariate_sumstats(Xtrain, Ytrain, standardize=standardize, 
                                             standardize.response=FALSE, mc.cores=nthreads)
grid <- autoselect.mixsd(univ_sumstats, mult=sqrt(2))^2

###Compute prior with only canonical and data-driven matrices
S0_can <- compute_canonical_covs(ncol(Ytrain), singletons=TRUE, hetgrid=c(0))
res_cor <- mashr::estimate_null_correlation_simple(mashr::mash_set_data(univ_sumstats$Bhat, univ_sumstats$Shat))
S0_datadriven <- compute_data_driven_covs(univ_sumstats, subset_thresh=0.05, n_pcs=3, 
                                          non_singleton=FALSE, Gamma=res_cor)
S0_full <- c(S0_can, S0_datadriven)
S0_full <- expand_covs(S0_full, grid, zeromat=TRUE)

###Fit grop-lasso to initialize mr.mash
if(nthreads>1){
  library(doMC)
  registerDoMC(nthreads)
  paral <- TRUE
} else {
  paral <- FALSE
}

cvfit_glmnet <- cv.glmnet(x=Xtrain, y=Ytrain, family="mgaussian", alpha=1, standardize=standardize, parallel=paral)
coeff_glmnet <- coef(cvfit_glmnet, s="lambda.min")
Bhat_glmnet <- matrix(as.numeric(NA), nrow=ncol(Xtrain), ncol=ncol(Ytrain))
ahat_glmnet <- rep(as.numeric(NA), ncol(Ytrain))
for(i in 1:length(coeff_glmnet)){
  Bhat_glmnet[, i] <- as.vector(coeff_glmnet[[i]])[-1]
  ahat_glmnet[i] <- as.vector(coeff_glmnet[[i]])[1]
}
Yhat_glmnet <- drop(predict(cvfit_glmnet, newx=Xtest, s="lambda.min"))
prop_nonzero_glmnet <- sum(Bhat_glmnet[, 1]!=0)/ncol(Xtrain)

###Fit mr.mash
w0 <- c((1-prop_nonzero_glmnet), rep(prop_nonzero_glmnet/(length(S0_full)-1), (length(S0_full)-1)))
fit_mrmash <- mr.mash(Xtrain, Ytrain, S0_full, w0=w0, update_w0=TRUE, update_w0_method="EM", tol=1e-2,
                      convergence_criterion="ELBO", compute_ELBO=TRUE, standardize=standardize, 
                      verbose=TRUE, update_V=TRUE, update_V_method="diagonal", e=1e-8,
                      mu1_init=Bhat_glmnet, w0_threshold=w0_threshold, nthreads=nthreads)
Yhat_mrmash <- predict(fit_mrmash, Xtest)


###Accuracy of prediction
cat("g-lasso mean scaled RMSE:\n")
print(mean(compute_accuracy(Ytest, Yhat_glmnet)$scaled_rmse))

cat("mr.mash mean scaled RMSE:\n")
print(mean(compute_accuracy(Ytest, Yhat_mrmash)$scaled_rmse))

cat("g-lasso individual tissue scaled RMSE:\n")
print(compute_accuracy(Ytest, Yhat_glmnet)$scaled_rmse)

cat("mr.mash individual tissue scaled RMSE:\n")
print(compute_accuracy(Ytest, Yhat_mrmash)$scaled_rmse)


