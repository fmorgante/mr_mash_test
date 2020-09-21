compute_accuracy <- function(Y, Yhat) {
  bias <- rep(as.numeric(NA), ncol(Y))
  r2 <- rep(as.numeric(NA), ncol(Y))
  rmse <- rep(as.numeric(NA), ncol(Y))
  
  for(i in 1:ncol(Y)){
    fit  <- lm(Y[, i] ~ Yhat[, i])
    bias[i] <- coef(fit)[2] 
    r2[i] <- summary(fit)$r.squared
    rmse[i] <- sqrt(mean((Y[, i] - Yhat[, i])^2))
  }
  
  return(list(bias=bias, r2=r2, rmse=rmse, scaled_rmse=rmse/matrixStats::colSds(Y)))
}
