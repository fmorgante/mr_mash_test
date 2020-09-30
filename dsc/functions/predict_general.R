predict.general <- function(B, intercept, newx){
  if(is.matrix(intercept))
    intercept <- drop(intercept)
  return(mr.mash.alpha:::addtocols(newx %*% B, intercept))
}