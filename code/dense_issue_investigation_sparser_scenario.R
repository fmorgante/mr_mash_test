library(mr.mash.alpha)
library(glmnet)

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

results <- vector("list", 5)

t <- 0

for(seed in c(1, 140, 303, 1024, 1999)){
  t <- t+1
  
  ###Set seed
  set.seed(seed)
  
  ###Set parameters
  n <- 1000
  p <- 1000
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
  test_set <- sort(sample(x=c(1:n), size=round(n*0.5), replace=FALSE))
  Ytrain <- out$Y[-test_set, ]
  Xtrain <- out$X[-test_set, ]
  Ytest <- out$Y[test_set, ]
  Xtest <- out$X[test_set, ]
  
  ###Build prior
  grid <- seq(0.1, 2.1, 0.2)
  S0 <- mr.mash.alpha:::compute_cov_canonical(ncol(Ytrain), singletons=TRUE, hetgrid=c(0, 0.25, 0.5, 0.75, 1), grid, zeromat=TRUE)
  
  ###Fit group-lasso
  cvfit_glmnet <- cv.glmnet(x=Xtrain, y=Ytrain, family="mgaussian", alpha=1)
  coeff_glmnet <- coef(cvfit_glmnet, s="lambda.min")
  Bhat_glmnet <- matrix(as.numeric(NA), nrow=p, ncol=r)
  ahat_glmnet <- vector("numeric", r)
  for(i in 1:length(coeff_glmnet)){
    Bhat_glmnet[, i] <- as.vector(coeff_glmnet[[i]])[-1]
    ahat_glmnet[i] <- as.vector(coeff_glmnet[[i]])[1]
  }
  Yhat_test_glmnet <- drop(predict(cvfit_glmnet, newx=Xtest, s="lambda.min"))
  
  ###Fit mr.mash updating free V and updating w0
  fit_mrmash <- mr.mash(Xtrain, Ytrain, S0, tol=1e-2, convergence_criterion="ELBO", update_w0=TRUE, 
                        update_w0_method="EM", compute_ELBO=TRUE, standardize=TRUE, verbose=FALSE, update_V=TRUE,
                        mu1_init = Bhat_glmnet, w0_threshold=1e-8)
  Yhat_test_mrmash <- predict(fit_mrmash, Xtest)
  
  ###Fit mr.mash fixed V and updating w0
  fit_mrmash_trueV <- mr.mash(Xtrain, Ytrain, S0, tol=1e-2, convergence_criterion="ELBO", update_w0=TRUE, 
                              update_w0_method="EM", compute_ELBO=TRUE, standardize=TRUE, verbose=FALSE, update_V=FALSE,
                              mu1_init = Bhat_glmnet, w0_threshold=1e-8, V=out$V)
  Yhat_test_mrmash_trueV <- predict(fit_mrmash_trueV, Xtest)
  
  ###Fit mr.mash with fixed V and w0 to the truth
  w0 <- rep(0, length(S0))
  w0[c(104, 111)] <- c(0.05, 0.95)
  fit_mrmash_trueV_truew0 <- mr.mash(Xtrain, Ytrain, S0, w0=w0, tol=1e-2, convergence_criterion="ELBO", update_w0=FALSE, 
                                     update_w0_method="EM", compute_ELBO=TRUE, standardize=TRUE, verbose=FALSE, update_V=FALSE,
                                     mu1_init=Bhat_glmnet, w0_threshold=1e-8, V=out$V)
  Yhat_test_mrmash_trueV_truew0 <- predict(fit_mrmash_trueV_truew0, Xtest)
  
  ###Prediction performance
  mse <- rbind(glmnet=compute_accuracy(Ytest, Yhat_test_glmnet)$mse,
               mrmash=compute_accuracy(Ytest, Yhat_test_mrmash)$mse,
               mrmash_trueV=compute_accuracy(Ytest, Yhat_test_mrmash_trueV)$mse,
               mrmash_trueV_truew0=compute_accuracy(Ytest, Yhat_test_mrmash_trueV_truew0)$mse)
  r2 <- rbind(glmnet=compute_accuracy(Ytest, Yhat_test_glmnet)$r2,
              mrmash=compute_accuracy(Ytest, Yhat_test_mrmash)$r2,
              mrmash_trueV=compute_accuracy(Ytest, Yhat_test_mrmash_trueV)$r2,
              mrmash_trueV_truew0=compute_accuracy(Ytest, Yhat_test_mrmash_trueV_truew0)$r2)
  bias <- rbind(glmnet=compute_accuracy(Ytest, Yhat_test_glmnet)$bias,
                mrmash=compute_accuracy(Ytest, Yhat_test_mrmash)$bias,
                mrmash_trueV=compute_accuracy(Ytest, Yhat_test_mrmash_trueV)$bias,
                mrmash_trueV_truew0=compute_accuracy(Ytest, Yhat_test_mrmash_trueV_truew0)$bias)
  
  ###ELBO
  elbo <- rbind(mrmash=fit_mrmash$ELBO, mrmash_trueV=fit_mrmash_trueV$ELBO, mrmash_trueV_truew0=fit_mrmash_trueV_truew0$ELBO)
  
  ###V
  V_true <- out$V
  V_mrmash <- fit_mrmash$V
  
  results[[t]] <- list(mse=mse, r2=r2, bias=bias, elbo=elbo, V_true=V_true, V_mrmash=V_mrmash)
}

saveRDS(results, "/project2/mstephens/fmorgante/mr_mash_test/output/test_sparse_issue.rds")
