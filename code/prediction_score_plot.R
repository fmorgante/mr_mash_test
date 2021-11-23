###Set options and load libraries
options(stringsAsFactors = FALSE)

library(optparse)
library(dplyr)

###Parse command line arguments
parser <- OptionParser()
parser <- add_option(parser, c("--score_dir"), type="character")
parser <- add_option(parser, c("--analysis_units"), type="character")
parser <- add_option(parser, c("--files_suffix"), type="character")
outparse <- parse_args(parser)

dir <- outparse$score_dir
suffix <- outparse$files_suffix
manifest <- outparse$analysis_units


# dir <- "/project/mstephens/fmorgante/mr_mash_test/output/gtex_mr_mash_analysis/prediction_score"
# suffix <- "GTEx_V8_score"
# manifest <- read.table("/project/mstephens/fmorgante/mr_mash_test/data/gtex-v8-manifest-2ormore-tissues-nopath-nosuffix-batch1.txt", header=FALSE)

###Create data frames to store the results
rmse_enet <- data.frame()
rmse_mrmash <- data.frame()

###Loop through genes in the manifest
for(i in 1:nrow(manifest)){
  ##Load scores data
  gene <- manifest[i,]
  dat <- readRDS(paste0(dir, "/", gene, ".", suffix, ".rds"))
  
  if(!is.null(dat)){
    ##Reshape scaled RMSE results
    rmse_enet_1 <- data.frame(matrix(dat$mean_scaled_rmse$enet, nrow=1))
    rmse_mrmash_1 <- data.frame(matrix(dat$mean_scaled_rmse$mrmash_first, nrow=1))
    colnames(rmse_enet_1) <- colnames(rmse_mrmash_1) <- dat$mean_scaled_rmse$tissue
    
    ##Add current results to the data frame
    rmse_enet <- bind_rows(rmse_enet, rmse_enet_1)
    rmse_mrmash <- bind_rows(rmse_mrmash, rmse_mrmash_1)
  } else {
    message(gene, " has an empty score file.")
  }
}

###Compute difference in RMSE between mr.mash and enet
rmse_diff <- as.matrix(rmse_mrmash) - as.matrix(rmse_enet)

###Plot results
##Pooled tissues
pdf("test_rmse_diff.pdf")
hist(rmse_diff)
dev.off()

##By tissue
rmse_diff_by_tissue_plot <- vector("list", ncol(rmse_diff))
for(j in 1:ncol(rmse_diff)){
  rmse_diff_by_tissue_plot[[j]] <- hist(rmse_diff[,j])
}

