# setwd("/project2/mstephens/fmorgante/software/mr.mash.alpha/inst/code")

# source("../../R/misc.R")
# source("../../R/mr_mash_simple.R")
library(mr.mash.alpha)
library(glmnet)
library(ggplot2)
library(cowplot)

###Set options
options(stringsAsFactors = FALSE)

###Function to compute accuracy
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

###Function to compute rmse (relative to group_lasso)
compute_rmse <- function(dat){
  dat <- transform(dat, experiment=paste(repl, resp, sep="-"))
  t <- 0
  for (i in unique(dat$experiment)) {
    t <- t+1
    rmse_data  <- dat[which(dat$experiment == i), ]
    mse_base <- rmse_data[which(rmse_data$method=="group_lasso"), "value"]
    rmse_data$value <- rmse_data$value/mse_base
    
    if(t>1){
      rmse_data_tot <- rbind(rmse_data_tot, rmse_data)
    } else if(t==1){
      rmse_data_tot <- rmse_data
    }
  }
  
  rmse_data_tot$experiment <- NULL
  
  return(rmse_data_tot)
}

###Loop over seeds
RNGversion("3.5.0")

###Set seed
set.seed(123)

###Set parameters
n <- 550
p <- 100
p_causal <- 50
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
test_set <- sort(sample(x=c(1:n), size=500, replace=FALSE))
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
Yhat_test_glmnet <- drop(predict(cvfit_glmnet, newx=Xtest, s="lambda.min"))

###Define prior
univ_sumstats <- mr.mash.alpha:::get_univariate_sumstats(Xtrain, Ytrain, standardize=FALSE, standardize.response=FALSE)
grid <- mr.mash.alpha:::autoselect.mixsd(c(univ_sumstats$Bhat), c(univ_sumstats$Shat), mult=sqrt(2))^2
S0 <- mr.mash.alpha:::compute_cov_canonical(ncol(Ytrain), singletons=TRUE, hetgrid=c(0, 0.25, 0.5, 0.75, 1), grid, zeromat=TRUE)

###Fit mr.mash fixing V to the truth
# w0 <- rep(1/length(S0), length(S0))
# fit_mrmash_simple_fixV <- mr_mash_simple(Xtrain, Ytrain, out$V, S0, w0, Bhat_glmnet, 200,
#                                          tol=1e-3, update_w0=TRUE, update_V=FALSE, verbose=TRUE)

fit_mrmash_fixV <- mr.mash.alpha:::mr.mash(Xtrain, Ytrain, S0=S0, tol=1e-3, convergence_criterion="mu1", update_w0=TRUE,
                      update_w0_method="EM", compute_ELBO=TRUE, standardize=FALSE, verbose=FALSE, update_V=FALSE,
                      mu1_init = Bhat_glmnet, w0_threshold=0, V=out$V)
Yhat_test_mrmash_fixV <- predict(fit_mrmash_fixV, Xtest)

###Fit mr.mash fixing V and g to the truth
closest_comp <- grid[which.min(abs(0.9-grid))]
S0fix <- list(S0_1=matrix(closest_comp, r, r), S0_2=matrix(0, r, r))
w0fix <- c(0.5, 0.5)

# fit_mrmash_simple_fixVandg <- mr_mash_simple(Xtrain, Ytrain, out$V, S0fix, w0fix, Bhat_glmnet, 200,
#                                          tol=1e-3, update_w0=FALSE, update_V=FALSE, verbose=TRUE)

fit_mrmash_fixVandg <- mr.mash(Xtrain, Ytrain, S0=S0fix, w0=w0fix, tol=1e-3, convergence_criterion="mu1", update_w0=FALSE,
                           update_w0_method="EM", compute_ELBO=TRUE, standardize=FALSE, verbose=FALSE, update_V=FALSE,
                           mu1_init = Bhat_glmnet, w0_threshold=0, V=out$V)
Yhat_test_mrmash_fixVandg <- predict(fit_mrmash_fixVandg, Xtest)

###Fit mr.mash fixing g to the truth
closest_comp <- grid[which.min(abs(0.9-grid))]
S0fix <- list(S0_1=matrix(closest_comp, r, r), S0_2=matrix(0, r, r))
w0fix <- c(0.5, 0.5)

# fit_mrmash_simple_g <- mr_mash_simple(Xtrain, Ytrain, cov(Ytrain - Xtrain%*%Bhat_glmnet), S0fix, w0fix, Bhat_glmnet, 200,
#                                       tol=1e-3, update_w0=FALSE, update_V=TRUE, verbose=TRUE)

fit_mrmash_fixg <- mr.mash(Xtrain, Ytrain, S0=S0fix, w0=w0fix, tol=1e-3, convergence_criterion="mu1", update_w0=FALSE,
                               update_w0_method="EM", compute_ELBO=TRUE, standardize=FALSE, verbose=FALSE, update_V=TRUE,
                               mu1_init = Bhat_glmnet, w0_threshold=0)
Yhat_test_mrmash_fixg <- predict(fit_mrmash_fixg, Xtest)

###Plot coeffcients
layout(matrix(c(1, 1, 2, 2,
                1, 1, 2, 2,
                3, 3, 4, 4,
                3, 3, 4, 4), 4, 4, byrow = TRUE))

print(plot(out$B, Bhat_glmnet, main="group-lasso", xlab="True coefficients", ylab="Estimated coefficients"))
print(abline(0, 1))
print(plot(out$B, fit_mrmash_fixV$mu1, main="mr.mash\nTrue V", xlab="True coefficients", ylab="Estimated coefficients"))
print(abline(0, 1))
print(plot(out$B, fit_mrmash_fixg$mu1, main="mr.mash\nTrue g", xlab="True coefficients", ylab="Estimated coefficients"))
print(abline(0, 1))
print(plot(out$B, fit_mrmash_fixVandg$mu1, main="mr.mash\nTrue V and g", xlab="True coefficients", ylab="Estimated coefficients"))
print(abline(0, 1))

###Compute mse of coefficients
compute_accuracy(out$B, Bhat_glmnet)$mse
compute_accuracy(out$B, fit_mrmash_fixV$mu1)$mse
compute_accuracy(out$B, fit_mrmash_fixVandg$mu1)$mse
compute_accuracy(out$B, fit_mrmash_fixg$mu1)$mse

###Compute prediction performance metrics
compute_accuracy(Ytest, Yhat_test_glmnet)$mse
compute_accuracy(Ytest, Yhat_test_mrmash_fixV)$mse
compute_accuracy(Ytest, Yhat_test_mrmash_fixVandg)$mse
compute_accuracy(Ytest, Yhat_test_mrmash_fixg)$mse

compute_accuracy(Ytest, Yhat_test_glmnet)$bias
compute_accuracy(Ytest, Yhat_test_mrmash_fixV)$bias
compute_accuracy(Ytest, Yhat_test_mrmash_fixVandg)$bias
compute_accuracy(Ytest, Yhat_test_mrmash_fixg)$bias

compute_accuracy(Ytest, Yhat_test_glmnet)$r2
compute_accuracy(Ytest, Yhat_test_mrmash_fixV)$r2
compute_accuracy(Ytest, Yhat_test_mrmash_fixVandg)$r2
compute_accuracy(Ytest, Yhat_test_mrmash_fixg)$r2
