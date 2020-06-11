library(mr.mash.alpha)
library(glmnet)
library(ggplot2)
library(cowplot)

setwd("/project2/mstephens/fmorgante/mr_mash_test")

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

###Create quantities to stores the results
value <- vector("numeric", 150)
method <- vector("character", 150)
repl <- vector("numeric", 150)
resp <- vector("numeric", 150)

bias <- data.frame(method, repl, resp, value)
r2 <- data.frame(method, repl, resp, value)
mse <- data.frame(method, repl, resp, value)
intercept <- data.frame(method, repl, resp, value)

z <- 0

###Loop over seeds
RNGversion("3.5.0")
for(seed in c(1, 65, 23, 5009, 2857, 33, 579, 319, 498, 1040)){
  ###Set seed
  set.seed(seed)

  ###Indices needed to place the results in the correct position  
  z <- z+1
  if(z==1){
    j <- 0
  } else {
    j <- (z-1)*25    
  }

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
  Yhat_test_glmnet <- drop(predict(cvfit_glmnet, newx=Xtest, s="lambda.min"))
  
  ###Define prior
  univ_sumstats <- mr.mash.alpha:::get_univariate_sumstats(Xtrain, Ytrain, standardize=FALSE, standardize.response=FALSE)
  grid <- mr.mash.alpha:::autoselect.mixsd(c(univ_sumstats$Bhat), c(univ_sumstats$Shat), mult=sqrt(2))^2
  S0 <- mr.mash.alpha:::compute_cov_canonical(ncol(Ytrain), singletons=TRUE, hetgrid=c(0, 0.25, 0.5, 0.75, 1), grid, zeromat=TRUE)
  
  ###Fit mr.mash estimating V 
  fit_mrmash <- mr.mash(Xtrain, Ytrain, S0, tol=1e-2, convergence_criterion="ELBO", update_w0=TRUE, 
                        update_w0_method="EM", compute_ELBO=TRUE, standardize=FALSE, verbose=FALSE, update_V=TRUE,
                        mu1_init = Bhat_glmnet, w0_threshold=1e-8)
  Yhat_test_mrmash <- predict(fit_mrmash, Xtest)
  
  ###Fit mr.mash fixing V to the truth
  fit_mrmash_fixV <- mr.mash(Xtrain, Ytrain, S0, tol=1e-2, convergence_criterion="ELBO", update_w0=TRUE, 
                             update_w0_method="EM", compute_ELBO=TRUE, standardize=FALSE, verbose=FALSE, update_V=FALSE,
                             mu1_init = Bhat_glmnet, w0_threshold=1e-8, V=out$V)
  Yhat_test_mrmash_fixV <- predict(fit_mrmash_fixV, Xtest)
  
  ###Fit mr.mash fixing V and g to the truth
  closest_comp <- grid[which.min(abs(0.9-grid))]
  S0fix <- list(S0_1=matrix(closest_comp, r, r), S0_2=matrix(0, r, r))
  w0fix <- c(0.5, 0.5)
  fit_mrmash_fixVandg <- mr.mash(Xtrain, Ytrain, S0=S0fix, w0=w0fix, tol=1e-2, convergence_criterion="ELBO", update_w0=FALSE, 
                                 update_w0_method="EM", compute_ELBO=TRUE, standardize=FALSE, verbose=FALSE, update_V=FALSE,
                                 mu1_init = Bhat_glmnet, w0_threshold=1e-8, V=out$V)
  Yhat_test_mrmash_fixVandg <- predict(fit_mrmash_fixVandg, Xtest)
  
  ###Fit mr.mash fixing g to the truth
  fit_mrmash_fixg <- mr.mash(Xtrain, Ytrain, S0=S0fix, w0=w0fix, tol=1e-2, convergence_criterion="ELBO", update_w0=FALSE, 
                             update_w0_method="EM", compute_ELBO=TRUE, standardize=FALSE, verbose=FALSE, update_V=TRUE,
                             mu1_init = Bhat_glmnet, w0_threshold=1e-8)
  Yhat_test_mrmash_fixg <- predict(fit_mrmash_fixg, Xtest)
  
  

  ###Bias (slope of Y~Yhat)
  bias[(j+1):(j+25), "value"] <- c(compute_accuracy(Ytest, Yhat_test_glmnet)$bias,
                        compute_accuracy(Ytest, Yhat_test_mrmash)$bias,
                        compute_accuracy(Ytest, Yhat_test_mrmash_fixV)$bias,
                        compute_accuracy(Ytest, Yhat_test_mrmash_fixVandg)$bias,
                        compute_accuracy(Ytest, Yhat_test_mrmash_fixg)$bias)
  bias[(j+1):(j+25), "method"] <- rep(c("group_lasso", "mr_mash", "mr_mash_trueV",
                                        "mr_mash_trueVandg", "mr_mash_trueg"), each=5)
  bias[(j+1):(j+25), "repl"] <- rep(z, 25)
  bias[(j+1):(j+25), "resp"] <- rep(1:5, times=5)
  
  ###r2
  r2[(j+1):(j+25), "value"] <- c(compute_accuracy(Ytest, Yhat_test_glmnet)$r2,
                             compute_accuracy(Ytest, Yhat_test_mrmash)$r2,
                             compute_accuracy(Ytest, Yhat_test_mrmash_fixV)$r2,
                             compute_accuracy(Ytest, Yhat_test_mrmash_fixVandg)$r2,
                             compute_accuracy(Ytest, Yhat_test_mrmash_fixg)$r2)
  r2[(j+1):(j+25), "method"] <- rep(c("group_lasso", "mr_mash", "mr_mash_trueV",
                                      "mr_mash_trueVandg", "mr_mash_trueg"), each=5)
  r2[(j+1):(j+25), "repl"] <- rep(z, 25)
  r2[(j+1):(j+25), "resp"] <- rep(1:5, times=5)
  
  ###MSE
  mse[(j+1):(j+25), "value"] <-  c(compute_accuracy(Ytest, Yhat_test_glmnet)$mse,
                               compute_accuracy(Ytest, Yhat_test_mrmash)$mse,
                               compute_accuracy(Ytest, Yhat_test_mrmash_fixV)$mse,
                               compute_accuracy(Ytest, Yhat_test_mrmash_fixVandg)$mse,
                               compute_accuracy(Ytest, Yhat_test_mrmash_fixg)$mse)
  mse[(j+1):(j+25), "method"] <- rep(c("group_lasso", "mr_mash", "mr_mash_trueV",
                                       "mr_mash_trueVandg", "mr_mash_trueg"), each=5)
  mse[(j+1):(j+25), "repl"] <- rep(z, 25)
  mse[(j+1):(j+25), "resp"] <- rep(1:5, times=5)
  
  
  ###Coefficients
  pdf(paste0("analysis/dense_issue_investigation2/coefficients_sim", z, ".pdf"))
  layout(matrix(c(1, 1, 2, 2,
                  1, 1, 2, 2,
                  3, 3, 4, 4,
                  3, 3, 4, 4), 4, 4, byrow = TRUE))

  print(plot(out$B, Bhat_glmnet, main="group-lasso", xlab="True coefficients", ylab="Estimated coefficients"))
  print(abline(0, 1))
  print(plot(out$B, fit_mrmash$mu1, main="mr.mash\nEstimated V", xlab="True coefficients", ylab="Estimated coefficients"))
  print(abline(0, 1))
  print(plot(out$B, fit_mrmash_fixV$mu1, main="mr.mash\nTrue V", xlab="True coefficients", ylab="Estimated coefficients"))
  print(abline(0, 1))
  print(plot(out$B, fit_mrmash_fixVandg$mu1, main="mr.mash\nTrue V and g", xlab="True coefficients", ylab="Estimated coefficients"))
  print(abline(0, 1))
  print(plot(out$B, fit_mrmash_fixg$mu1, main="mr.mash\nTrue g", xlab="True coefficients", ylab="Estimated coefficients"))
  print(abline(0, 1))
  
  dev.off()
  
  ###Intercept
  intercept[(j+1):(j+25), "value"] <- c(ahat_glmnet, fit_mrmash$intercept, fit_mrmash_fixV$intercept, 
                                        fit_mrmash_fixVandg$intercept, fit_mrmash_fixg$intercept)
  intercept[(j+1):(j+25), "method"] <- rep(c("group_lasso", "mr_mash", "mr_mash_trueV",
                                             "mr_mash_trueVandg", "mr_mash_trueg"), each=5)
  intercept[(j+1):(j+25), "repl"] <- rep(z, 25)
  intercept[(j+1):(j+25), "resp"] <- rep(1:5, times=5)
  
  
  ###ELBO
  # rbind(mr_mash=fit_mrmash$ELBO, mr_mash_trueV=fit_mrmash_fixV$ELBO)
  
  ###Mixture weights
  # sort(fit_mrmash$w0, decreasing=TRUE)
  # sort(fit_mrmash_fixV$w0, decreasing=TRUE)
  
  ###Residual covariance
  # out$V
  # fit_mrmash$V
}

###Compute relative (to group_lasso) MSE
rmse <- compute_rmse(mse)

###Add factor version of resp
rmse$resp_fac <- as.factor(rmse$resp)
bias$resp_fac <- as.factor(bias$resp)
r2$resp_fac <- as.factor(r2$resp)
intercept$resp_fac <- as.factor(intercept$resp)

###Plots
p_rmse <- ggplot(rmse, aes_string(x = "resp_fac", y = "value", fill = "method")) +
  geom_boxplot(color = "black", outlier.size = 1, width = 0.85) +
  labs(x = "Response", y = "Error", title = "MSE", fill="Method") +
  theme_cowplot(font_size = 18) +
  theme(plot.title = element_text(hjust = 0.5))

pdf("analysis/dense_issue_investigation2/rmse.pdf")
p_rmse
dev.off()

p_bias <- ggplot(bias, aes_string(x = "resp_fac", y = "value", fill = "method")) +
  geom_boxplot(color = "black", outlier.size = 1, width = 0.85) +
  labs(x = "Response", y = "Slope", title = "Bias", fill="Method") +
  geom_hline(yintercept=1, linetype="dashed", size=1) +
  theme_cowplot(font_size = 18) +
  theme(plot.title = element_text(hjust = 0.5))

pdf("analysis/dense_issue_investigation2/bias.pdf")
p_bias
dev.off()

p_r2 <- ggplot(r2, aes_string(x = "resp_fac", y = "value", fill = "method")) +
  geom_boxplot(color = "black", outlier.size = 1, width = 0.85) +
  labs(x = "Response", y = "Accuracy", title = "r2", fill="Method") +
  theme_cowplot(font_size = 18) +
  theme(plot.title = element_text(hjust = 0.5))

pdf("analysis/dense_issue_investigation2/r2.pdf")
p_r2
dev.off()

p_intercept <- ggplot(intercept, aes_string(x = "resp_fac", y = "value", fill = "method")) +
  geom_boxplot(color = "black", outlier.size = 1, width = 0.85) +
  labs(x = "Response", y = "Intercept", title = "Intercept", fill="Method") +
  geom_hline(yintercept=1, linetype="dashed", size=1) +
  theme_cowplot(font_size = 18) +
  theme(plot.title = element_text(hjust = 0.5))

pdf("analysis/dense_issue_investigation2/intercept.pdf")
p_intercept
dev.off()


