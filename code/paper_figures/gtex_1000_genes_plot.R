options(stringsAsFactors = FALSE)

###Load libraries
library(reshape2)
library(dplyr)
library(matrixStats)
library(ggplot2)
library(cowplot)

dir <- "../../output/gtex_mr_mash_analysis/prediction_score"
suffix <- "GTEx_V8_score"
manifest <- read.table("../../data/gtex-v8-manifest-2ormore-tissues-nopath-nosuffix-batch1.txt", header=FALSE)
sample_size <- read.csv("../../data/gtex-v8-sample-size-by-tissue.csv")

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

###Compute relative improvement of mr.mash over the enet based on RMSE
rmse_diff <- (as.matrix(rmse_mrmash) - as.matrix(rmse_enet))/as.matrix(rmse_enet)
medians <- colMedians(rmse_diff, na.rm=TRUE)
means <- data.frame(tissue=colnames(rmse_diff), mean_RMSE_diff=colMeans2(rmse_diff, na.rm=TRUE))
rmse_diff <- rmse_diff[, order(medians, decreasing = FALSE)]
rmse_diff_melt <- reshape2::melt(rmse_diff, value.name="RMSE_diff")
rmse_diff_melt <- rmse_diff_melt[, -1]
colnames(rmse_diff_melt)[1] <- "Tissue"


###Plot of relative improvement of mr.mash over the enet based on RMSE across tissues
par(mar = c(13.1,5.4,4,2) + 0.1)
p <- ggplot(rmse_diff_melt, aes_string(x = "Tissue", y = "RMSE_diff")) +
  geom_boxplot(color = "black", fill="gold", outlier.size = 1, width = 0.85) +
  stat_summary(fun=mean, geom="point", shape=18, size=2, color="black", fill="yellow") +
  labs(x = "Tissue", y = "Difference in RMSE relative to e-net") +
  geom_hline(yintercept=0, linetype="dotted", size=1, color = "red") +
  theme_cowplot(font_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust=1, size=8))

ggsave("../../analysis/paper_figures/gtex_1000_genes_mrmash_vs_enet_rmse.pdf", plot=p, device="pdf", units="in", height=8, width=11)


###Compute mean relative improvement of mr.mash over the enet based on RMSE by tissue
rmse_diff_means <- data.frame(Tissue=colnames(rmse_diff), mean_RMSE_diff=colMeans2(rmse_diff, na.rm=TRUE))
rmse_diff_means_sample_size <- merge(rmse_diff_means, sample_size, by="Tissue", all.x=TRUE)

###Plot of relative improvement of mr.mash over the enet based on RMSE vs sample size by tissue
p_size <- ggplot(rmse_diff_means_sample_size, aes(x=X..RNASeq.and.Genotyped.samples, y=mean_RMSE_diff)) + 
                   geom_point(shape=16) +
                   geom_smooth(method='lm', formula= y~x) +
                   labs(x = "Sample size", y = "Difference in RMSE relative to e-net") + 
                   theme_cowplot(font_size = 16)

ggsave("../../analysis/paper_figures/gtex_1000_genes_mrmash_vs_enet_rmse_vs_sample_size.pdf", plot=p_size, device="pdf", units="in", height=8, width=8)



