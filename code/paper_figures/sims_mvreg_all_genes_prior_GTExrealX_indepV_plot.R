###Set options
options(stringsAsFactors=FALSE)

###Load libraries
library(ggplot2)
library(cowplot)
library(scales)

dsc_plots <- readRDS("sims_paper_figures_inter/")

###Set some quantities used in the following plots
colors <- c("cyan", "skyblue", "dodgerblue", "mediumblue", "limegreen", "green", "gold", "orange", "red", "firebrick", "darkmagenta", "mediumpurple")

###mr.mash -- different priors
##Create plot
methods_chosen <- c("mr_mash_can", "mr_mash_data", "mr_mash_both")
metric_chosen <- "scaled_rrmse"

p_priors <- ggplot(dsc_plots[which(dsc_plots$score_metric==metric_chosen & dsc_plots$method_fac %in% methods_chosen), ], 
                  aes_string(x = "response_fac", y = "score_value", fill = "method_fac")) +
  geom_boxplot(color = "black", outlier.size = 1, width = 0.85) +
  scale_fill_manual(values = colors[7:9], labels = c("canonical", "data-driven", "both")) +
  labs(x = "Tissue", y = "RMSE relative to both-drop", fill="Method") +
  geom_hline(yintercept=1, linetype="dotted", size=1) +
  theme_cowplot(font_size = 20) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("test.pdf", plot=p_priors, device="pdf", units="in", height=10, width=15)

###mr.mash vs other methods
##Create plot
methods_chosen <- c("mlasso", "mtlasso", "enet")
metric_chosen <- "scaled_rrmse"

p_methods <- ggplot(dsc_plots[which(dsc_plots$score_metric==metric_chosen & dsc_plots$method_fac %in% methods_chosen), ], 
                  aes_string(x = "response_fac", y = "score_value", fill = "method_fac")) +
  geom_boxplot(color = "black", outlier.size = 1, width = 0.85) +
  scale_fill_manual(values = colors[c(3,7,12)], labels = c("g-lasso", "smt-lasso", "e-net")) +
  labs(x = "Tissue", y = expression(paste("RMSE relative to ", italic("mr.mash"))), fill="Method") +
  geom_hline(yintercept=1, linetype="dotted", size=1) +
  theme_cowplot(font_size = 20) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("test1.pdf", plot=p_methods, device="pdf", units="in", height=10, width=15)


