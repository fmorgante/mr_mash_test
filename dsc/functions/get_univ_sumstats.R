get_univ_sumstats <- function(X, Y, standardize, zscores, nthreads){
  p <- ncol(X)
  r <- ncol(Y)
  
  univ_sumstats <- mr.mash.alpha::compute_univariate_sumstats(X=X, Y=Y, standardize=standardize, 
                                                              standardize.response=FALSE, 
                                                              mc.cores=nthreads)
  rownames(univ_sumstats$Bhat) <- paste0("X", 1:p)
  rownames(univ_sumstats$Shat) <- paste0("X", 1:p)
  colnames(univ_sumstats$Bhat) <- paste0("Y", 1:r)
  colnames(univ_sumstats$Shat) <- paste0("Y", 1:r)
  
  out <- list(bhat=univ_sumstats$Bhat, shat=univ_sumstats$Shat)
  
  if(zscores){
    Zscores <- univ_sumstats$Bhat/univ_sumstats$Shat
    rownames(Zscores) <- paste0("X", 1:p)
    colnames(Zscores) <- paste0("Y", 1:r)
    
    out$zhat <- Zscores
  }
  
  return(out)
}

