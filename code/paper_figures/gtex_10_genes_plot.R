options(stringsAsFactors = FALSE)

###Load libraries
library(dplyr)
library(ggplot2)
library(cowplot)

dir <- "../../output/gtex_mr_mash_analysis/prediction_score_mtlasso"
suffix <- "GTEx_V8_score"
manifest <- read.table("../../data/gtex-v8-manifest-2ormore-tissues-nopath-nosuffix-mtlasso.txt", header=FALSE)

###Create data frames to store the results
rmse_mtlasso <- data.frame()
rmse_mrmash <- data.frame()

###Loop through genes in the manifest
for(i in 1:nrow(manifest)){
  ##Load scores data
  gene <- manifest[i,]
  dat <- readRDS(paste0(dir, "/", gene, ".", suffix, ".rds"))
  
  if(!is.null(dat)){
    ##Reshape scaled RMSE results
    rmse_mtlasso_1 <- data.frame(matrix(dat$mean_scaled_rmse$mtlasso, nrow=1))
    rmse_mrmash_1 <- data.frame(matrix(dat$mean_scaled_rmse$mrmash_first, nrow=1))
    colnames(rmse_mtlasso_1) <- colnames(rmse_mrmash_1) <- dat$mean_scaled_rmse$tissue
    rownames(rmse_mtlasso_1) <- rownames(rmse_mrmash_1) <- gene
    
    ##Add current results to the data frame
    rmse_mtlasso <- bind_rows(rmse_mtlasso, rmse_mtlasso_1)
    rmse_mrmash <- bind_rows(rmse_mrmash, rmse_mrmash_1)
  } else {
    message(gene, " has an empty score file.")
  }
}

p <- vector("list", nrow(rmse_mtlasso))

###Loop through genes that did not fail and make the plot objects
for(i in 1:nrow(rmse_mtlasso)){
  min_lim <- min(c(unlist(rmse_mrmash[i, ]), unlist(rmse_mtlasso[i, ])))
  max_lim <- max(c(unlist(rmse_mrmash[i, ]), unlist(rmse_mtlasso[i, ])))
  dat_plot <- data.frame(mrmash_first=unlist(rmse_mrmash[i, ]), mtlasso=unlist(rmse_mtlasso[i, ]))
  
  p[[i]] <- ggplot(dat_plot, aes(x=mtlasso, y=mrmash_first)) + geom_point(shape=1) + xlim(min_lim, max_lim) + ylim(min_lim, max_lim) +
                   labs(x = "smt-lasso RMSE", y = expression(paste(italic("mr.mash"), " RMSE")), title=rownames(rmse_mrmash)[i]) + 
                   geom_abline(intercept = 0, slope = 1) + 
                   theme_cowplot(font_size = 12) +
                   theme(plot.title = element_text(hjust = 0.5, size=12), axis.text = element_text(size=8))  
}

###Arrange the plots into a grid
p_grid <- plot_grid(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]],
                    p[[6]], p[[7]], p[[8]], p[[9]], NULL,
                    p[[10]], NULL,
                    ncol = 3)

ggsave("../../analysis/paper_figures/gtex_10_genes_mrmash_vs_mtlasso_rmse.pdf", plot=p_grid, device="pdf", units="in", height=11, width=11)
