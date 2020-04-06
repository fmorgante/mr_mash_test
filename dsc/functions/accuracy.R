accuracy <- function(Y, Yhat) {
  bias <- rep(NA, ncol(Y))
  r2 <- rep(NA, ncol(Y))
  mse <- rep(NA, ncol(Y))
  
  for(i in 1:ncol(Y)){
    fit  <- lm(Y[, i] ~ Yhat[, i])
    bias[i] <- coef(fit)[2] 
    r2[i] <- summary(fit)$r.squared
    mse[i] <- mean((Y[, i] - Yhat[, i])^2)
  }
  
  return(list(bias=bias, r2=r2, mse=mse))
}