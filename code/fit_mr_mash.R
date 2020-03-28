########################
## Fabio Morgante
## 3/27/2020
## Test mr.mash for prediction
########################

###Load libraries
library(mr.mash.alpha)
library(optparse)

###Set options
options(stringsAsFactors = FALSE)

###Set seed
set.seed(1)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--n"), type="integer")
parser <- add_option(parser, c("--p"), type="integer")
parser <- add_option(parser, c("--p_causal"), type="integer")
parser <- add_option(parser, c("--r"), type="integer")
parser <- add_option(parser, c("--pve"), type="numeric")
parser <- add_option(parser, c("--sigma_offdiag"), type="numeric")
parser <- add_option(parser, c("--sigma_scale"), type="numeric")
parser <- add_option(parser, c("--gamma_offdiag"), type="numeric")
parser <- add_option(parser, c("--gamma_scale"), type="numeric")
parser <- add_option(parser, c("--V_offdiag"), type="numeric")
parser <- add_option(parser, c("--V_scale"), type="numeric")
parser <- add_option(parser, c("--prop_testset"), type="numeric", default=0.2)
parser <- add_option(parser, c("--update_w0"), type="logical", default=TRUE)
parser <- add_option(parser, c("--update_w0_method"), type="character", default="EM")
parser <- add_option(parser, c("--standardize"), type="logical", default=TRUE)
parser <- add_option(parser, c("--verbose"), type="logical", default=TRUE)
parser <- add_option(parser, c("--update_V"), type="logical", default=TRUE)
parser <- add_option(parser, c("--outdir"), type="character", default="./")
outparse <- parse_args(parser)

###Simulate V, B, X and Y
dat <- mr.mash.alpha:::simulate_mr_mash_data(outparse$n, outparse$p, outparse$p_causal, outparse$r, intercepts=rep(1, outparse$r),
                                             pve=outparse$pve, Sigma_cor_offdiag=outparse$sigma_offdiag, Sigma_scale=outparse$sigma_scale,
                                             Gamma_cor_offdiag=outparse$gamma_offdiag, Gamma_scale=outparse$gamma_scale,
                                             V_cor_offdiag=outparse$V_offdiag, V_offdiag_scale=outparse$V_scale)

###Split data in training and test sets
test_set <- sort(sample(x=c(1:outparse$n), size=round(outparse$n*outparse$prop_testset), replace=FALSE))
Ytrain <- dat$Y[-test_set, ]
Xtrain <- dat$X[-test_set, ]
Ytest <- dat$Y[test_set, ]
Xtest <- dat$X[test_set, ]

###Build the mixture prior
grid <- seq(0.1, 2, 0.2)
S0 <- mr.mash.alpha:::compute_cov_canonical(outparse$r, singletons=TRUE, hetgrid=c(0, 0.25, 0.5, 0.75, 0.99), grid, zeromat=TRUE)

###Fit model in the training set
fit <- mr.mash(Xtrain, Ytrain, S0, tol=1e-8, update_w0=outparse$update_w0, update_w0_method=outparse$update_w0_method,
               compute_ELBO=TRUE, standardize=outparse$standardize, verbose=outparse$verbose, update_V=outparse$update_V,
               version="Rcpp", e=1e-8)

###Predict phenotypes in the test set 
Yhat_test <- predict(fit, Xtest)

###Prepare and save the output
output <- list(params=outparse, inputs=dat, test_set=test_set, Ytrain=Ytrain, Xtrain=Xtrain, Ytest=Ytest, Xtest=Xtest, fit=fit,
               Yhat_test=Yhat_test)
if(!outparse$update_w0){
  update_w0_method <- "NA"
}else{
  update_w0_method <- outparse$update_w0_method
}
saveRDS(output, paste0(outparse$outdir, "fit_mr_mash_", "n", outparse$n, "_p", outparse$p, "_p_caus", outparse$p_causal, "_r", outparse$r,
                       "_pve", outparse$pve, "_sigmaoffdiag", outparse$sigma_offdiag, "_sigmascale", outparse$sigma_scale,
                       "_gammaoffdiag", outparse$gamma_offdiag, "_gammascale", outparse$gamma_scale,
                       "_Voffdiag", outparse$V_offdiag, "_Vscale", outparse$V_scale, "_updatew0", outparse$update_w0,
                       "_updatew0", outparse$update_w0, "_updatew0method", update_w0_method, "_updateV", outparse$update_V, ".rds"))

###Print session info
print(sessionInfo())
