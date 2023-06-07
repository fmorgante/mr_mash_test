options(stringsAsFactors=FALSE)

library(dplyr)
library(stringr)


###List all the mr.mash output files
files <- list.files(path="/project/mstephens/fmorgante/mr_mash_test/output/gtex_mr_mash_analysis/prediction/fold_1", 
                    pattern = "\\.first_pass.rds$", full.names=TRUE)

i = 0

###Load results for each gene and compile them in a single data.frame
for (f in files) {
  i = i+1    
  
  dat = readRDS(f)$model$w0
  if (is.null(dat)) {
    message(paste("Dataset", f, "has no valid w0 quantity"))
    next
  }
  if(i > 1){
    weights <- bind_rows(weights, dat)
  } else {
    weights <- dat
  }
}

###Compute mean weights across genes
weights1 <- as.data.frame(weights)
weights1[is.na(weights1)] <- 0 #assign NA to mixture components that were dropped during model fit
weights1_means <- colMeans(weights1)

###Read in estimated data-driven matrices
mats <- readRDS("/project/mstephens/fmorgante/mr_mash_test/output/gtex_mr_mash_analysis/data_driven_matrices/output/fold_1.ted_unconstrained.rds")

mats_names <- names(mats$U)
mats_names <- paste0(mats_names, sep="_")

###Create data.frame to store results
comp <- vector(mode="character", length(mats_names))
perc_nonnull_weight <- vector(mode="numeric", length(mats_names))
res <- data.frame(comp, perc_nonnull_weight)

j <- 0

###Loop over matrices and compute the percentage of non-null weight
nonnull_weight <- 1 - weights1_means["null"]

for(m in mats_names){
  j <- j+1
  tokeep <- stringr::str_detect(names(weights1_means), m)
  res[j, 1] <- m
  res[j, 2] <- (sum(weights1_means[tokeep])/nonnull_weight)*100
}

###Save the results
saveRDS(res, "/project/mstephens/fmorgante/mr_mash_test/perc_nonnull_mixture_weights_across_genes_fold_1.rds")
