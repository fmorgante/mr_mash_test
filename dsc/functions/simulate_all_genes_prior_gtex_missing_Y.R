simulate_data_all_genes_prior_gtex_missing_Y <- function(X, p_causal, r, r_causal, pve, B_cor, B_scale, 
                                                          w, V_cor, testset_index, prop_miss){
  
  if(!is.numeric(p_causal)){
    p_causal <- sample(x=eval(parse(text=p_causal)), size=1)
  }
  
  dat <- simulate_mr_mash_data_from_given_X(X=X, p_causal=p_causal, r=r, r_causal=r_causal, intercepts=rep(1, r),
                                            pve=pve, B_cor=B_cor, B_scale=B_scale, w=w, V_cor=V_cor)
  
  colnames(dat$Y) <- paste0("Y", seq(1, r))
  rownames(dat$Y) <- rownames(X)
  colnames(dat$X) <- colnames(X)
  rownames(dat$X) <- rownames(X)
  
  test_set <- readRDS(testset_index)
  Ytrain <- dat$Y[!(rownames(dat$Y) %in% test_set), ]
  Xtrain <- dat$X[!(rownames(dat$X) %in% test_set), ]
  Ytest <- dat$Y[test_set, ]
  Xtest <- dat$X[test_set, ]
  
  Y_miss <- assign_NAs_random(Ytrain, prop_miss)
  all_miss <- which(apply(Y_miss, 1, function(x, r){sum(is.na(x))==r}, r))
  if(length(all_miss)>0){
    Xtrain <- Xtrain[-all_miss, ]
    Ytrain <- Y_miss[-all_miss, ]
  } else {
    Ytrain <- Y_miss
  }
  
  return(list(Xtrain=Xtrain, Ytrain=Ytrain, Xtest=Xtest, Ytest=Ytest, B_true=dat$B))
}