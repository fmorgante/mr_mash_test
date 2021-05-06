compute_coeff_mse <- function(B, Bhat) {
  mse <- rep(as.numeric(NA), ncol(B))
  
    for(i in 1:ncol(B)){
      mse[i] <- mean((B[, i] - Bhat[, i])^2)
    }

  return(mse)
}
