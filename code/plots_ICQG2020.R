################
###          ###
### ACCURACY ###
###          ###
################

###Set options
options(stringsAsFactors=FALSE)

###Load libraries
library(dscrutils)
library(ggplot2)
library(cowplot)
library(scales)


#############################
### Pooled across tissues ###
#############################

###Function to convert dscquery output from list to data.frame suitable for plotting
convert_dsc_to_dataframe_gtex <- function(dsc){
  ###Data.frame to store the results after convertion
  dsc_df <- data.frame()
  
  ###Get length of list elements 
  n_elem <- length(dsc$DSC)
  
  ###Loop through the dsc list
  for(i in 1:n_elem){
    ##Prepare vectors making up the final data frame
    r_scalar <- dsc$simulate.r[i]
    response <- 1:r_scalar
    pve <- rep(dsc$simulate.pve[i], times=r_scalar)
    simulate <- rep(dsc$simulate[i], times=r_scalar)
    fit <- rep(dsc$fit[i], times=r_scalar)
    score <- rep(dsc$score[i], times=r_scalar)
    score.err <- dsc$score.err[[i]]
    timing <- rep(dsc$fit.time[i], times=r_scalar)
    
    ##Build the data frame
    df <- data.frame(response=response, pve=pve, scenario=simulate, method=fit, 
                     score_metric=score, score_value=score.err, time=timing)
    dsc_df <- rbind(dsc_df, df)
  }
  
  ###Get number of genes
  n_methods <- length(unique(dsc$fit))
  n_metrics <- length(unique(dsc$score))
  n_genes <- n_elem/n_methods/n_metrics
  
  ###Compute replicates
  reptot <- rep(rep(1:n_genes, each=dsc$simulate.r[1]), n_methods*n_metrics)
  
  dsc_df <- data.frame(rep=reptot, dsc_df)
  
  return(dsc_df)
}

###Function to compute rrmse 
compute_rrmse <- function(dsc_plot, base_level="mr_mash_em_can", log10_scale=FALSE){
  dsc_plot <- transform(dsc_plot, experiment=paste(rep, response, scenario, sep="-"))
  t <- 0
  for (i in unique(dsc_plot$experiment)) {
    t <- t+1
    rmse_data  <- dsc_plot[which(dsc_plot$experiment == i & dsc_plot$score_metric=="scaled_rmse"), ]
    mse_mr_mash_consec_em <- rmse_data[which(rmse_data$method==base_level), "score_value"]
    if(!log10_scale)
      rmse_data$score_value <- rmse_data$score_value/mse_mr_mash_consec_em
    else
      rmse_data$score_value <- log10(rmse_data$score_value/mse_mr_mash_consec_em)
    rmse_data$score_metric <- "scaled_rrmse"
    if(t>1){
      rmse_data_tot <- rbind(rmse_data_tot, rmse_data)
    } else if(t==1){
      rmse_data_tot <- rmse_data
    }
  }
  
  rmse_data_tot$experiment <- NULL
  
  return(rmse_data_tot)
}

###Set some quantities used in the following plots
colors <- c("cyan", "skyblue", "limegreen", "green", "gold")

###Scenario
path_files <- c("../output/mvreg_all_genes_prior_GTExrealX_indepV_indepB_r10",
                "../output/mvreg_all_genes_prior_GTExrealX_indepV_sharedB_r10",
                "../output/mvreg_all_genes_prior_GTExrealX_indepV_indepB_r10_lowPVE",
                "../output/mvreg_all_genes_prior_GTExrealX_indepV_sharedB_r10_lowPVE",
                "../output/mvreg_all_genes_prior_GTExrealX_indepV_B_1causalrespr10",
                "../output/mvreg_all_genes_prior_GTExrealX_indepV_sharedB_3causalrespr10")

for(path_file in path_files){
  
  ###Load the dsc results
  dsc_out <- dscquery(path_file, 
                      c("simulate.r", "simulate.pve", "simulate", "fit", "score", "score.err", "fit.time"), 
                      groups="fit: mr_mash_em_can, mr_mash_em_data, mr_mash_em_dataAndcan, mr_mash_em_dataAndcan_dropcomp, mlasso, mtlasso, enet",
                      verbose=FALSE)
  
  ###Convert from list to data.frame for plotting
  dsc_plots <- convert_dsc_to_dataframe_gtex(dsc_out)
  
  ###Remove unused methods
  dsc_plots <- dsc_plots[which(dsc_plots$method!="mr_mash_em_can" & dsc_plots$method!="mr_mash_em_dataAndcan" & 
                                 dsc_plots$method!="mr_mash_em_dataAndcan_dropcomp"), ]
  
  ###Compute rmse score (relative to mr_mash_consec_em) and add it to the data
  rrmse_dat <- compute_rrmse(dsc_plots, base_level="mr_mash_em_data")
  dsc_plots <- rbind(dsc_plots, rrmse_dat)
  
  ###Remove undesired scores and keep only methods wanted
  dsc_plots <- dsc_plots[which(dsc_plots$score_metric!="scaled_rmse"), ]
  dsc_plots <- dsc_plots[which(dsc_plots$method!="mr_mash_em_data"), ]
  
  ###Create factor version of method
  dsc_plots$method_fac <- factor(dsc_plots$method, levels=c("mlasso", "mtlasso", "enet"),
                                 labels=c("g-lasso", "smt-lasso", "e-net"))
  
  ###Create plots
  p_rrmse <- ggplot(dsc_plots, aes_string(x = "method_fac", y = "score_value", fill = "method_fac")) +
    #geom_boxplot(color = "black", outlier.size = 1, width = 0.85) +
    Ipaper::geom_boxplot2(color = "black", width.errorbar = 0, width = 0.85) +
    scale_fill_manual(values = colors, guide=FALSE) +
    labs(x = "", y = "RMSE relative to mr.mash") +
    geom_hline(yintercept=1, linetype="dotted", size=1) +
    theme_cowplot(font_size = 16) +
    theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1))
  
  ###Plot
  #print(p_rrmse)
  simscen_name <- unlist(strsplit(path_file, "/"))[3]
  if(path_file=="../output/mvreg_all_genes_prior_GTExrealX_indepV_B_1causalrespr10" | 
     path_file=="../output/mvreg_all_genes_prior_GTExrealX_indepV_sharedB_3causalrespr10"){
    simscen_name_plot <- paste0("../analysis/", simscen_name, "_pooled.pdf")
  }else{
    simscen_name_plot <- paste0("../analysis/", simscen_name, ".pdf")
  }
  ggsave(simscen_name_plot, p_rrmse, width = 8.5, height = 5.5, units = "in")
}


rm(list=ls())


##############################
### Each tissue separately ###
##############################

###Function to convert dscquery output from list to data.frame suitable for plotting
convert_dsc_to_dataframe_gtex <- function(dsc){
  ###Data.frame to store the results after convertion
  dsc_df <- data.frame()
  
  ###Get length of list elements 
  n_elem <- length(dsc$DSC)
  
  ###Loop through the dsc list
  for(i in 1:n_elem){
    ##Prepare vectors making up the final data frame
    r_scalar <- dsc$simulate.r[i]
    response <- 1:r_scalar
    pve <- rep(dsc$simulate.pve[i], times=r_scalar)
    simulate <- rep(dsc$simulate[i], times=r_scalar)
    fit <- rep(dsc$fit[i], times=r_scalar)
    score <- rep(dsc$score[i], times=r_scalar)
    score.err <- dsc$score.err[[i]]
    timing <- rep(dsc$fit.time[i], times=r_scalar)
    
    ##Build the data frame
    df <- data.frame(response=response, pve=pve, scenario=simulate, method=fit, 
                     score_metric=score, score_value=score.err, time=timing)
    dsc_df <- rbind(dsc_df, df)
  }
  
  ###Get number of genes
  n_methods <- length(unique(dsc$fit))
  n_metrics <- length(unique(dsc$score))
  n_genes <- n_elem/n_methods/n_metrics
  
  ###Compute replicates
  reptot <- rep(rep(1:n_genes, each=dsc$simulate.r[1]), n_methods*n_metrics)
  
  dsc_df <- data.frame(rep=reptot, dsc_df)
  
  return(dsc_df)
}

###Function to compute rrmse 
compute_rrmse <- function(dsc_plot, base_level="mr_mash_em_can", log10_scale=FALSE){
  dsc_plot <- transform(dsc_plot, experiment=paste(rep, response, scenario, sep="-"))
  t <- 0
  for (i in unique(dsc_plot$experiment)) {
    t <- t+1
    rmse_data  <- dsc_plot[which(dsc_plot$experiment == i & dsc_plot$score_metric=="scaled_rmse"), ]
    mse_mr_mash_consec_em <- rmse_data[which(rmse_data$method==base_level), "score_value"]
    if(!log10_scale)
      rmse_data$score_value <- rmse_data$score_value/mse_mr_mash_consec_em
    else
      rmse_data$score_value <- log10(rmse_data$score_value/mse_mr_mash_consec_em)
    rmse_data$score_metric <- "scaled_rrmse"
    if(t>1){
      rmse_data_tot <- rbind(rmse_data_tot, rmse_data)
    } else if(t==1){
      rmse_data_tot <- rmse_data
    }
  }
  
  rmse_data_tot$experiment <- NULL
  
  return(rmse_data_tot)
}

###Set some quantities used in the following plots
colors <- c("cyan", "skyblue", "limegreen", "green", "gold")

###List of files to use
path_files <- c("../output/mvreg_all_genes_prior_GTExrealX_indepV_sharedB_2blocksr10",
                "../output/mvreg_all_genes_prior_GTExrealX_indepV_B_1causalrespr10",
                "../output/mvreg_all_genes_prior_GTExrealX_indepV_sharedB_3causalrespr10",
                "../output/mvreg_all_genes_prior_GTExrealX_indepV_sharedB_2blocksr10_lowPVE",
                "../output/mvreg_all_genes_prior_GTExrealX_indepV_sharedB_2blocksr10_diffPVE",
                "../output/mvreg_all_genes_prior_GTExrealX_indepV_B_1causalrespr10_lowPVE",
                "../output/mvreg_all_genes_prior_GTExrealX_indepV_sharedB_3causalrespr10_lowPVE")


for(path_file in path_files){
  
  ###Load the dsc results
  dsc_out <- dscquery(path_file, 
                      c("simulate.r", "simulate.pve", "simulate", "fit", "score", "score.err", "fit.time"), 
                      groups="fit: mr_mash_em_can, mr_mash_em_data, mr_mash_em_dataAndcan, mr_mash_em_dataAndcan_dropcomp, mlasso, mtlasso, enet",
                      verbose=FALSE)
  
  ###Convert from list to data.frame for plotting
  dsc_plots <- convert_dsc_to_dataframe_gtex(dsc_out)
  
  ###Remove unused methods
  dsc_plots <- dsc_plots[which(dsc_plots$method!="mr_mash_em_can" & dsc_plots$method!="mr_mash_em_dataAndcan" & 
                                 dsc_plots$method!="mr_mash_em_dataAndcan_dropcomp"), ]
  
  ###Compute rmse score (relative to mr_mash_consec_em) and add it to the data
  rrmse_dat <- compute_rrmse(dsc_plots, base_level="mr_mash_em_data")
  dsc_plots <- rbind(dsc_plots, rrmse_dat)
  
  ###Remove undesired scores and keep only methods wanted
  dsc_plots <- dsc_plots[which(dsc_plots$score_metric!="scaled_rmse"), ]
  dsc_plots <- dsc_plots[which(dsc_plots$score_metric!="r2"), ]
  dsc_plots <- dsc_plots[which(dsc_plots$method!="mr_mash_em_data"), ]
  
  ###Create factor version of method
  dsc_plots$method_fac <- factor(dsc_plots$method, levels=c("mlasso", "mtlasso", "enet"),
                                 labels=c("g-lasso", "smt-lasso", "e-net"))
  

  ###Create factor version of response
  dsc_plots$response_fac <- as.factor(dsc_plots$response)

  ###Create plots
  p_rrmse <- ggplot(dsc_plots, aes_string(x = "response_fac", y = "score_value", fill = "method_fac")) +
    #geom_boxplot(color = "black", outlier.size = 1, width = 0.85) +
    Ipaper::geom_boxplot2(color = "black", width.errorbar = 0, width = 0.85) +
    scale_fill_manual(values = colors) + #, guide=guide_legend(nrow=2,byrow=TRUE)
    labs(x = "Tissue", y = "RMSE relative to mr.mash", fill="Method") +
    geom_hline(yintercept=1, linetype="dotted", size=1) +
    theme_cowplot(font_size = 16) +
    theme(legend.position="top", legend.justification='center') 

  ###Plot
  #print(p_rrmse)
  simscen_name <- unlist(strsplit(path_file, "/"))[3]
  simscen_name_plot <- paste0("../analysis/", simscen_name, ".pdf")
  ggsave(simscen_name_plot, p_rrmse, width = 8.5, height = 5.5, units = "in")
}


rm(list=ls())


###############
###         ###
### RUNTIME ###
###         ###
###############

###Function to convert dscquery output from list to data.frame suitable for plotting
convert_dsc_to_dataframe_gtex <- function(dsc){
  ###Data.frame to store the results after convertion
  dsc_df <- data.frame()

  ###Get length of list elements
  n_elem <- length(dsc$DSC)

  ###Loop through the dsc list
  for(i in 1:n_elem){
    ##Prepare vectors making up the final data frame
    r_scalar <- dsc$simulate.r[i]
    response <- 1:r_scalar
    pve <- rep(dsc$simulate.pve[i], times=r_scalar)
    simulate <- rep(dsc$simulate[i], times=r_scalar)
    fit <- rep(dsc$fit[i], times=r_scalar)
    score <- rep(dsc$score[i], times=r_scalar)
    score.err <- dsc$score.err[[i]]
    timing <- rep(dsc$fit.time[i], times=r_scalar)

    ##Build the data frame
    df <- data.frame(response=response, pve=pve, scenario=simulate, method=fit,
                     score_metric=score, score_value=score.err, time=timing)
    dsc_df <- rbind(dsc_df, df)
  }

  ###Get number of genes
  n_methods <- length(unique(dsc$fit))
  n_metrics <- length(unique(dsc$score))
  n_genes <- n_elem/n_methods/n_metrics

  ###Compute replicates
  reptot <- rep(rep(1:n_genes, each=dsc$simulate.r[1]), n_methods*n_metrics)

  dsc_df <- data.frame(rep=reptot, dsc_df)

  return(dsc_df)
}

###List of files to use
path_files <- c("../output/mvreg_all_genes_prior_GTExrealX_indepV_sharedB_r10",
                "../output/mvreg_all_genes_prior_GTExrealX_indepV_sharedB_2blocksr10_diffPVE",
                "../output/mvreg_all_genes_prior_GTExrealX_indepV_B_1causalrespr10",
                "../output/mvreg_all_genes_prior_GTExrealX_indepV_sharedB_3causalrespr10")

###Create data frame to store the results
dsc_plots_time <- data.frame()

###Loop over files
for(path_file in path_files){
  ###Load the dsc results
  dsc_out <- dscquery(path_file,
                      c("simulate.r", "simulate.pve", "simulate", "fit", "score", "score.err", "fit.time"),
                      groups="fit: mr_mash_em_can, mr_mash_em_data, mr_mash_em_dataAndcan, mr_mash_em_dataAndcan_dropcomp, mlasso, mtlasso, enet",
                      verbose=FALSE)

  ###Convert from list to data.frame for plotting
  dsc_plots <- convert_dsc_to_dataframe_gtex(dsc_out)
  
  ###Remove unused methods
  dsc_plots <- dsc_plots[which(dsc_plots$method!="mr_mash_em_can" & dsc_plots$method!="mr_mash_em_dataAndcan" & 
                                 dsc_plots$method!="mr_mash_em_dataAndcan_dropcomp"), ]
  
  ###Create time dataset
  dsc_plots_time <- rbind(dsc_plots_time, dsc_plots[which(dsc_plots$response==1 & dsc_plots$score_metric=="scaled_rmse"),
                                                    -which(colnames(dsc_plots) %in% c("score_metric", "score_value", "response"))])
}

###Set some quantities used in the following plots
colors_time <- c("gold", "cyan", "skyblue", "limegreen", "green")

###Create factor version of method
dsc_plots_time$method_fac <- factor(dsc_plots_time$method, levels=c("mr_mash_em_data", "mlasso", "mtlasso", "enet"),
                               labels=c("mr.mash", "g-lasso", "smt-lasso", "e-net"))

###Create plots
p_time <- ggplot(dsc_plots_time, aes_string(x = "method_fac", y = "time", fill = "method_fac")) +
  #geom_boxplot(color = "black", outlier.size = 1, width = 0.85) +
  Ipaper::geom_boxplot2(color = "black", width.errorbar = 0, width = 0.85) +
  scale_fill_manual(values = colors_time, guide=FALSE) +
  scale_y_continuous(trans="log2", breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x))) +
  labs(x = "", y = "Elapsed time (seconds) in log2 scale") +
  theme_cowplot(font_size = 16) +
  theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1))

###Plot
#print(p_time)
ggsave("../analysis/pooled_time.pdf", p_time, width = 8.5, height = 5.5, units = "in")


rm(list=ls())


##########
### r2 ###
##########

###Function to convert dscquery output from list to data.frame suitable for plotting
convert_dsc_to_dataframe_gtex <- function(dsc){
  ###Data.frame to store the results after convertion
  dsc_df <- data.frame()
  
  ###Get length of list elements 
  n_elem <- length(dsc$DSC)
  
  ###Loop through the dsc list
  for(i in 1:n_elem){
    ##Prepare vectors making up the final data frame
    r_scalar <- dsc$simulate.r[i]
    response <- 1:r_scalar
    pve <- rep(dsc$simulate.pve[i], times=r_scalar)
    simulate <- rep(dsc$simulate[i], times=r_scalar)
    fit <- rep(dsc$fit[i], times=r_scalar)
    score <- rep(dsc$score[i], times=r_scalar)
    score.err <- dsc$score.err[[i]]
    timing <- rep(dsc$fit.time[i], times=r_scalar)
    
    ##Build the data frame
    df <- data.frame(response=response, pve=pve, scenario=simulate, method=fit, 
                     score_metric=score, score_value=score.err, time=timing)
    dsc_df <- rbind(dsc_df, df)
  }
  
  ###Get number of genes
  n_methods <- length(unique(dsc$fit))
  n_metrics <- length(unique(dsc$score))
  n_genes <- n_elem/n_methods/n_metrics
  
  ###Compute replicates
  reptot <- rep(rep(1:n_genes, each=dsc$simulate.r[1]), n_methods*n_metrics)
  
  dsc_df <- data.frame(rep=reptot, dsc_df)
  
  return(dsc_df)
}

###Function to compute rrmse 
compute_rrmse <- function(dsc_plot, base_level="mr_mash_em_can", log10_scale=FALSE){
  dsc_plot <- transform(dsc_plot, experiment=paste(rep, response, scenario, sep="-"))
  t <- 0
  for (i in unique(dsc_plot$experiment)) {
    t <- t+1
    rmse_data  <- dsc_plot[which(dsc_plot$experiment == i & dsc_plot$score_metric=="scaled_rmse"), ]
    mse_mr_mash_consec_em <- rmse_data[which(rmse_data$method==base_level), "score_value"]
    if(!log10_scale)
      rmse_data$score_value <- rmse_data$score_value/mse_mr_mash_consec_em
    else
      rmse_data$score_value <- log10(rmse_data$score_value/mse_mr_mash_consec_em)
    rmse_data$score_metric <- "scaled_rrmse"
    if(t>1){
      rmse_data_tot <- rbind(rmse_data_tot, rmse_data)
    } else if(t==1){
      rmse_data_tot <- rmse_data
    }
  }
  
  rmse_data_tot$experiment <- NULL
  
  return(rmse_data_tot)
}

###Set some quantities used in the following plots
colors <- c("gold", "cyan", "skyblue", "limegreen", "green")

###Load the dsc results
dsc_out <- dscquery("../output/mvreg_all_genes_prior_GTExrealX_indepV_sharedB_2blocksr10_diffPVE", 
                    c("simulate.r", "simulate.pve", "simulate", "fit", "score", "score.err", "fit.time"), 
                    groups="fit: mr_mash_em_can, mr_mash_em_data, mr_mash_em_dataAndcan, mr_mash_em_dataAndcan_dropcomp, mlasso, mtlasso, enet",
                    verbose=FALSE)
  
###Convert from list to data.frame for plotting
dsc_plots <- convert_dsc_to_dataframe_gtex(dsc_out)
  
###Remove unused methods
dsc_plots <- dsc_plots[which(dsc_plots$method!="mr_mash_em_can" & dsc_plots$method!="mr_mash_em_dataAndcan" & 
                               dsc_plots$method!="mr_mash_em_dataAndcan_dropcomp"), ]
  
###Remove undesired scores
dsc_plots <- dsc_plots[which(dsc_plots$score_metric=="r2"), ]
  
###Create factor version of method
dsc_plots$method_fac <- factor(dsc_plots$method, levels=c("mr_mash_em_data", "mlasso", "mtlasso", "enet"),
                               labels=c("mr.mash", "g-lasso", "smt-lasso", "e-net"))
  
###Create factor version of response
dsc_plots$response_fac <- as.factor(dsc_plots$response)
  
###Create plots
p_r2 <- ggplot(dsc_plots, aes_string(x = "response_fac", y = "score_value", fill = "method_fac")) +
  #geom_boxplot(color = "black", outlier.size = 1, width = 0.85) +
  Ipaper::geom_boxplot2(color = "black", width.errorbar = 0, width = 0.85) +
  scale_fill_manual(values = colors) + #, guide=guide_legend(nrow=2,byrow=TRUE)
  labs(x = "Tissue", y = expression(r^2), fill="Method") +
  geom_hline(yintercept=0.2, linetype="dotted", size=1) +
  geom_hline(yintercept=0.05, linetype="dashed", size=1) +
  theme_cowplot(font_size = 16) +
  theme(legend.position="top", legend.justification='center') 
  
###Plot
#print(p_r2)
ggsave("../analysis/mvreg_all_genes_prior_GTExrealX_indepV_sharedB_2blocksr10_diffPVE_r2.pdf", p_r2, width = 8.5, height = 5.5, units = "in")

