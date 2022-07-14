options(stringsAsFactors = FALSE)

###Load libraries
library(foreach)
library(doMC)
library(optparse)
library(matrixStats)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--n_cores"), type="integer")
outparse <- parse_args(parser)


###Functions to compute MAF and missing genotype rate
compute_maf <- function(geno){
  f <- mean(geno,na.rm = TRUE)/2
  return(min(f, 1-f))
}

compute_missing <- function(geno){
  miss <- sum(is.na(geno))/length(geno)
  return(miss)
}

mean_impute <- function(geno){
  f <- apply(geno, 2, function(x) mean(x,na.rm = TRUE))
  for (i in 1:length(f)) geno[,i][which(is.na(geno[,i]))] <- f[i]
  return(geno)
}

### Filter X matrix
filter_X <- function(X, missing_rate_thresh, maf_thresh, var_thresh) {
  rm_col <- which(apply(X, 2, compute_missing) > missing_rate_thresh)
  if (length(rm_col)) X <- X[, -rm_col]
  rm_col <- which(apply(X, 2, compute_maf) < maf_thresh)
  if (length(rm_col)) X <- X[, -rm_col]
  X <- mean_impute(X)
  #rm_col <- which(apply(X, 2, is_zero_variance))
  rm_col <- which(matrixStats::colVars(X) < var_thresh)
  if (length(rm_col)) X <- X[, -rm_col]
  return(X)
}

###Register cores
registerDoMC(outparse$n_cores)

###Load complete manifest
manifest <- read.table("../data/gtex-v8-manifest-full-X.txt", header=FALSE)

###Obtain indexes of genes with full X (i.e., 838 individuals)
res <- foreach(i = 1:nrow(manifest), .combine='rbind') %dopar% {
  X <- readRDS(manifest[i, 1])$X
  if(nrow(X)==838){
    X_filt <- filter_X(X, 0.05, 0.05, 0.05)
    data.frame(gene=paste(unlist(strsplit(unlist(strsplit(manifest[i, 1], "/"))[4], "[.]"))[1:2], collapse="."), 
               n_snps_total=ncol(X), n_snps_filtered=ncol(X_filt))
  }
}

###Write out results
write.table(res, "../data/gtex-v8-full-X-number-of-snps.txt", sep="\t", row.names = FALSE, col.names=TRUE, quote=FALSE)
