###Set options
options(stringsAsFactors=FALSE)

###Load libraries
library(ggplot2)
library(cowplot)
library(scales)

dsc_plots_1causalresp <- readRDS("../../output/sims_paper_figures_inter/mvreg_all_genes_prior_GTExrealX_indepV_B_1causalrespr10_indepReps_missing_Y_for_plotting.rds")
dsc_plots_indep <- readRDS("../../output/sims_paper_figures_inter/mvreg_all_genes_prior_GTExrealX_indepV_indepB_r10_indepReps_missing_Y_for_plotting.rds")
dsc_plots_shared <- readRDS("../../output/sims_paper_figures_inter/mvreg_all_genes_prior_GTExrealX_indepV_sharedB_r10_indepReps_missing_Y_for_plotting.rds")
dsc_plots_2blocks <- readRDS("../../output/sims_paper_figures_inter/mvreg_all_genes_prior_GTExrealX_indepV_sharedB_2blocksr10_diffPVE_indepReps_missing_Y_for_plotting.rds")
dsc_plots_3causalresp <- readRDS("../../output/sims_paper_figures_inter/mvreg_all_genes_prior_GTExrealX_indepV_sharedB_3causalrespr10_indepReps_missing_Y_for_plotting.rds")

###Set some quantities used in the following plots
colors <- c("cyan", "skyblue", "dodgerblue", "mediumblue", "limegreen", "green", "gold", "orange", "red", "firebrick", "darkmagenta", "mediumpurple")

####mr.mash vs other methods
###Error
##Create plot
methods_chosen <- c("mr_mash_em_data_mean_impute_enet", "mtlasso", "enet")
metric_chosen <- "scaled_rrmse"

#Filter data
dsc_plots_methods_1causalresp_filt <- dsc_plots_1causalresp[which(dsc_plots_1causalresp$score_metric==metric_chosen & dsc_plots_1causalresp$method %in% methods_chosen), ]
dsc_plots_methods_indep_filt <- dsc_plots_indep[which(dsc_plots_indep$score_metric==metric_chosen & dsc_plots_indep$method %in% methods_chosen), ]
dsc_plots_methods_shared_filt <- dsc_plots_shared[which(dsc_plots_shared$score_metric==metric_chosen & dsc_plots_shared$method %in% methods_chosen), ]  
dsc_plots_methods_2blocks_filt <- dsc_plots_2blocks[which(dsc_plots_2blocks$score_metric==metric_chosen & dsc_plots_2blocks$method %in% methods_chosen), ]
dsc_plots_methods_3causalresp_filt <- dsc_plots_3causalresp[which(dsc_plots_3causalresp$score_metric==metric_chosen & dsc_plots_3causalresp$method %in% methods_chosen), ]

#Extract y-axis limits
# ymin <- min(min(dsc_plots_1causalresp_filt$score_value), min(dsc_plots_indep_filt$score_value),
# 			min(dsc_plots_shared_filt$score_value), min(dsc_plots_2blocks_filt$score_value),
# 			min(dsc_plots_3causalresp_filt$score_value))
# 			
# ymax <- max(max(dsc_plots_1causalresp_filt$score_value), max(dsc_plots_indep_filt$score_value),
# 			max(dsc_plots_shared_filt$score_value), max(dsc_plots_2blocks_filt$score_value),
# 			max(dsc_plots_3causalresp_filt$score_value))

#Make plot objects for each scenario
p_methods_1causalresp <- ggplot(dsc_plots_methods_1causalresp_filt, 
                  aes_string(x = "response_fac", y = "score_value", fill = "method_fac")) +
  geom_boxplot(color = "black", outlier.size = 1, width = 0.85) +
#  ylim(ymin, ymax) +
  scale_fill_manual(values = colors[c(5,7,12)], labels = c(expression(paste(italic("mr.mash"), "(mean impute)")), "smt-lasso", "e-net")) +
  labs(x = "Tissue", y = expression(paste("RMSE relative to ", italic("mr.mash"))), title = "Mostly null", fill="Method") +
  geom_hline(yintercept=1, linetype="dotted", size=1, color = "red") +
  theme_cowplot(font_size = 16) +
  theme(plot.title = element_text(hjust = 0.5, size=14), legend.position="none")
  
p_methods_indep <- ggplot(dsc_plots_methods_indep_filt, 
                  aes_string(x = "response_fac", y = "score_value", fill = "method_fac")) +
  geom_boxplot(color = "black", outlier.size = 1, width = 0.85) +
#  ylim(ymin, ymax) +
  scale_fill_manual(values = colors[c(5,7,12)], labels = c(expression(paste(italic("mr.mash"), "(mean impute)")), "smt-lasso", "e-net")) +
  labs(x = "Tissue", y = expression(paste("RMSE relative to ", italic("mr.mash"))), title = "Independent effects", fill="Method") +
  geom_hline(yintercept=1, linetype="dotted", size=1, color = "red") +
  theme_cowplot(font_size = 16) +
  theme(plot.title = element_text(hjust = 0.5, size=14), legend.position="none")

p_methods_shared <- ggplot(dsc_plots_methods_shared_filt, 
                  aes_string(x = "response_fac", y = "score_value", fill = "method_fac")) +
  geom_boxplot(color = "black", outlier.size = 1, width = 0.85) +
#  ylim(ymin, ymax) +
  scale_fill_manual(values = colors[c(5,7,12)], labels = c(expression(paste(italic("mr.mash"), "(mean impute)")), "smt-lasso", "e-net")) +
  labs(x = "Tissue", y = expression(paste("RMSE relative to ", italic("mr.mash"))), title = "Equal effects", fill="Method") +
  geom_hline(yintercept=1, linetype="dotted", size=1, color = "red") +
  theme_cowplot(font_size = 16) +
  theme(plot.title = element_text(hjust = 0.5, size=14), legend.position="none")
  
p_methods_2blocks_inter <- ggplot(dsc_plots_methods_2blocks_filt, 
                  aes_string(x = "response_fac", y = "score_value", fill = "method_fac")) +
  geom_boxplot(color = "black", outlier.size = 1, width = 0.85) +
#  ylim(ymin, ymax) +
  scale_fill_manual(values = colors[c(5,7,12)], labels = c(expression(paste(italic("mr.mash"), "(mean impute)")), "smt-lasso", "e-net")) +
  labs(x = "Tissue", y = expression(paste("RMSE relative to ", italic("mr.mash"))), title = "Shared effects in subgroups", fill="Method") +
  geom_hline(yintercept=1, linetype="dotted", size=1, color = "red") +
  theme_cowplot(font_size = 16) +
  theme(plot.title = element_text(hjust = 0.5, size=14),
        legend.text.align = 0)
  
p_methods_2blocks <- p_methods_2blocks_inter + theme(legend.position="none")

p_methods_3causalresp <- ggplot(dsc_plots_methods_3causalresp_filt, 
                  aes_string(x = "response_fac", y = "score_value", fill = "method_fac")) +
  geom_boxplot(color = "black", outlier.size = 1, width = 0.85) +
#  ylim(ymin, ymax) +
  scale_fill_manual(values = colors[c(5,7,12)], labels = c(expression(paste(italic("mr.mash"), "(mean impute)")), "smt-lasso", "e-net")) +
  labs(x = "Tissue", y = expression(paste("RMSE relative to ", italic("mr.mash"))), title = "Equal effects + null", fill="Method") +
  geom_hline(yintercept=1, linetype="dotted", size=1, color = "red") +
  theme_cowplot(font_size = 16) +
  theme(plot.title = element_text(hjust = 0.5, size=14), legend.position="none")
  
#Extract legend
leg_methods <- get_legend(p_methods_2blocks_inter)

#Make the multi panel plot 
p_methods_all <- plot_grid(p_methods_shared, p_methods_indep, p_methods_2blocks,
			   p_methods_1causalresp, p_methods_3causalresp, 
			   leg_methods, labels = c('A', 'B', 'C', 'D', 'E'))

ggsave("../../analysis/paper_figures/mvreg_all_genes_prior_GTExrealX_indepV_missing_Y_mrmash_vs_others_rmse.pdf", plot=p_methods_all, device="pdf", units="in", height=8, width=11)


###Time
methods_chosen <- c("mr_mash_em_data_enet", "mr_mash_em_data_mean_impute_enet", "mtlasso", "enet")
metric_chosen <- "scaled_rrmse"

#Filter data
dsc_plots_1causalresp_methods_time <- dsc_plots_1causalresp[which(dsc_plots_1causalresp$response==1 & dsc_plots_1causalresp$score_metric==metric_chosen &
                                                          dsc_plots_1causalresp$method %in% methods_chosen), ]
 
dsc_plots_indep_methods_time <- dsc_plots_indep[which(dsc_plots_indep$response==1 & dsc_plots_indep$score_metric==metric_chosen &
                                                          dsc_plots_indep$method %in% methods_chosen), ]

dsc_plots_shared_methods_time <- dsc_plots_shared[which(dsc_plots_shared$response==1 & dsc_plots_shared$score_metric==metric_chosen &
                                                          dsc_plots_shared$method %in% methods_chosen), ]

dsc_plots_2blocks_methods_time <- dsc_plots_2blocks[which(dsc_plots_2blocks$response==1 & dsc_plots_2blocks$score_metric==metric_chosen &
                                                          dsc_plots_2blocks$method %in% methods_chosen), ]

dsc_plots_3causalresp_methods_time <- dsc_plots_3causalresp[which(dsc_plots_3causalresp$response==1 & dsc_plots_3causalresp$score_metric==metric_chosen &
                                                          dsc_plots_3causalresp$method %in% methods_chosen), ]

#Extract y-axis limits
ymin <- min(min(dsc_plots_1causalresp_methods_time$time), min(dsc_plots_indep_methods_time$time),
                      min(dsc_plots_shared_methods_time$time), min(dsc_plots_2blocks_methods_time$time),
                      min(dsc_plots_3causalresp_methods_time$time))

ymax <- max(max(dsc_plots_1causalresp_methods_time$time), max(dsc_plots_indep_methods_time$time),
                      max(dsc_plots_shared_methods_time$time), max(dsc_plots_2blocks_methods_time$time),
                      max(dsc_plots_3causalresp_methods_time$time))
                                       
p_methods_time_1causalresp <- ggplot(dsc_plots_1causalresp_methods_time, aes_string(x = "method_fac", y = "time", fill = "method_fac")) +
  geom_boxplot(color = "black", outlier.size = 1, width = 0.85) +
  scale_fill_manual(values = colors[c(9,5,7,12)], labels = c(expression(italic("mr.mash")), expression(paste(italic("mr.mash"), "(mean impute)")), "smt-lasso", "e-net")) +
  scale_y_continuous(limits=c(ymin,ymax), trans="log10", breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  labs(x = "", y = "Elapsed time (seconds)", title = "Mostly null", fill="Method") +
  theme_cowplot(font_size = 16) +
  theme(axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position="none", plot.title = element_text(hjust = 0.5, size=14))
        
p_methods_time_indep <- ggplot(dsc_plots_indep_methods_time, aes_string(x = "method_fac", y = "time", fill = "method_fac")) +
  geom_boxplot(color = "black", outlier.size = 1, width = 0.85) +
  scale_fill_manual(values = colors[c(9,5,7,12)], labels = c(expression(italic("mr.mash")), expression(paste(italic("mr.mash"), "(mean impute)")), "smt-lasso", "e-net")) +
  scale_y_continuous(limits=c(ymin,ymax), trans="log10", breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  labs(x = "", y = "Elapsed time (seconds)", title = "Independent effects", fill="Method") +
  theme_cowplot(font_size = 16) +
  theme(axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position="none", plot.title = element_text(hjust = 0.5, size=14))

p_methods_time_shared <- ggplot(dsc_plots_shared_methods_time, aes_string(x = "method_fac", y = "time", fill = "method_fac")) +
  geom_boxplot(color = "black", outlier.size = 1, width = 0.85) +
  scale_fill_manual(values = colors[c(9,5,7,12)], labels = c(expression(italic("mr.mash")), expression(paste(italic("mr.mash"), "(mean impute)")), "smt-lasso", "e-net")) +
  scale_y_continuous(limits=c(ymin,ymax), trans="log10", breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  labs(x = "", y = "Elapsed time (seconds)", title = "Equal effects", fill="Method") +  
  theme_cowplot(font_size = 16) +
  theme(axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position="none", plot.title = element_text(hjust = 0.5, size=14))

p_methods_time_2blocks_inter <- ggplot(dsc_plots_2blocks_methods_time, aes_string(x = "method_fac", y = "time", fill = "method_fac")) +
  geom_boxplot(color = "black", outlier.size = 1, width = 0.85) +
  scale_fill_manual(values = colors[c(9,5,7,12)], labels = c(expression(italic("mr.mash")), expression(paste(italic("mr.mash"), "(mean impute)")), "smt-lasso", "e-net")) +
  scale_y_continuous(limits=c(ymin,ymax), trans="log10", breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  labs(x = "", y = "Elapsed time (seconds)", title = "Shared effects in subgroups", fill="Method") +
  theme_cowplot(font_size = 16) +
  theme(axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.text.align = 0)
        
p_methods_time_2blocks <- p_methods_time_2blocks_inter + theme(legend.position="none", plot.title = element_text(hjust = 0.5, size=14))

p_methods_time_3causalresp <- ggplot(dsc_plots_3causalresp_methods_time, aes_string(x = "method_fac", y = "time", fill = "method_fac")) +
  geom_boxplot(color = "black", outlier.size = 1, width = 0.85) +
  scale_fill_manual(values = colors[c(9,5,7,12)], labels = c(expression(italic("mr.mash")), expression(paste(italic("mr.mash"), "(mean impute)")), "smt-lasso", "e-net")) +
  scale_y_continuous(limits=c(ymin,ymax), trans="log10", breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  labs(x = "", y = "Elapsed time (seconds)", title = "Equal effects + null", fill="Method") +
  theme_cowplot(font_size = 16) +
  theme(axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position="none", plot.title = element_text(hjust = 0.5, size=14))


#Extract legend
leg_methods_time <- get_legend(p_methods_time_2blocks_inter)

#Make the multi panel plot 
p_methods_time_all <- plot_grid(p_methods_time_shared, p_methods_time_indep, p_methods_time_2blocks,
				p_methods_time_1causalresp, p_methods_time_3causalresp, 
		  		leg_methods_time, labels = c('A', 'B', 'C', 'D', 'E'))

ggsave("../../analysis/paper_figures/mvreg_all_genes_prior_GTExrealX_indepV_missing_Y_mrmash_vs_others_runtime.pdf", plot=p_methods_time_all, device="pdf", units="in", height=8, width=11)
