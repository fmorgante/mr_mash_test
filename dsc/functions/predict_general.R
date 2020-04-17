predict.general <- function(B, intercept, newx){
  return(mr.mash.alpha:::addtocols(newx %*% B, intercept))
}