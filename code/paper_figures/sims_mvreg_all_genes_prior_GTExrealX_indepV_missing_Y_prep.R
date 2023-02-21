###Set options
options(stringsAsFactors=FALSE)

###Load libraries
library(dscrutils)
library(optparse)

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
    p_causal <- rep(dsc$simulate.p_causal[i], times=r_scalar)
    r <- rep(dsc$simulate.r[i], times=r_scalar)
    response <- 1:r_scalar
    pve <- rep(dsc$simulate.pve[i], times=r_scalar)
    simulate <- rep(dsc$simulate[i], times=r_scalar)
    fit <- rep(dsc$fit[i], times=r_scalar)
    score <- rep(dsc$score[i], times=r_scalar)
    score.err <- dsc$score.err[[i]]
    timing <- rep(dsc$fit.time[i], times=r_scalar)
    
    ##Build the data frame
    df <- data.frame(p_num_caus=p_causal, r=r, response=response, pve=pve,  
                     scenario=simulate, method=fit, score_metric=score, score_value=score.err, time=timing)
    dsc_df <- rbind(dsc_df, df)
  }
  
  ###Get number of genes
  n_methods <- length(unique(dsc$fit))
  n_metrics <- length(unique(dsc$score))
  n_genes <- n_elem/n_methods/n_metrics
  
  ###Compute replicates
  # beginning <- (repli*n_genes)-(n_genes-1)
  # ending <- repli*n_genes
  beginning <- 1
  ending <- n_genes
  reptot <- rep(rep(beginning:ending, each=dsc$simulate.r[1]), n_methods*n_metrics)
  
  dsc_df <- data.frame(subrep=reptot, dsc_df)
  
  return(dsc_df)
}

###Function to compute rmse (relative to a baseline)
compute_rrmse <- function(dsc_plot, baseline){
  dsc_plot <- transform(dsc_plot, experiment=paste(rep, subrep, response, scenario, sep="-"))
  t <- 0
  for (i in unique(dsc_plot$experiment)) {
    t <- t+1
    rmse_data  <- dsc_plot[which(dsc_plot$experiment == i & dsc_plot$score_metric=="scaled_rmse"), ]
    mse_mr_mash_consec_em <- rmse_data[which(rmse_data$method==baseline), "score_value"]
    rmse_data$score_value <- rmse_data$score_value/mse_mr_mash_consec_em
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

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--input"), type="character")
outparse <- parse_args(parser)
input <- outparse$input

###Loop over replicates
for(repl in 1:20){
  ###Load the dsc results
  dsc_out <- dscquery(paste0("../../output/", input, "/rep", repl), 
                       c("simulate.p_causal", "simulate.r", "simulate.r_causal", "simulate.pve", "simulate.B_cor", "simulate.w", 
                        "simulate.B_scale", "simulate.V_cor", "simulate", "fit", "score", "score.err", "fit.time"), 
                      groups="fit: mr_mash_em_data_mean_impute_enet, mr_mash_em_data_enet, mr_mash_em_can_enet, mr_mash_em_can_mlasso, mr_mash_em_data_mlasso, mr_mash_em_dataAndcan_mlasso, mr_mash_em_dataAndcan_dropcomp_mlasso, mtlasso, enet",
		      verbose=FALSE)
 
  ###Convert from list to data.frame for plotting
  if(repl==1){
    dsc_plots <- convert_dsc_to_dataframe_gtex(dsc_out)
    dsc_plots <- cbind(rep=repl, dsc_plots)
  } else if(repl>1){
    tmp <- convert_dsc_to_dataframe_gtex(dsc_out)
    tmp <- cbind(rep=repl, tmp)
    dsc_plots <- rbind(dsc_plots, tmp)
  }
}

###Compute rmse score and add it to the data
rrmse_dat <- compute_rrmse(dsc_plots, baseline="mr_mash_em_data_enet")
dsc_plots <- rbind(dsc_plots, rrmse_dat)

###Remove mse from scores and keep only methods wanted
#dsc_plots <- dsc_plots[which(dsc_plots$score_metric!="scaled_rmse" ), ]

###Create factor version of method
dsc_plots$method_fac <- factor(dsc_plots$method, levels=c("mr_mash_em_can_mlasso", "mr_mash_em_data_mlasso",
                                                          "mr_mash_em_dataAndcan_mlasso", "mr_mash_em_dataAndcan_dropcomp_mlasso", 
                                                          "mr_mash_em_can_enet", "mr_mash_em_data_enet", "mr_mash_em_data_mean_impute_enet", "mtlasso", "enet"),
                               labels=c("mr_mash_can_mlasso", "mr_mash_data_mlasso", "mr_mash_both_mlasso", 
                                        "mr_mash_both_drop_mlasso", "mr_mash_can_enet", "mr_mash_data_enet", 
                                        "mr_mash_data_mean_imp_enet", "mtlasso", "enet"))

###Create factor version of response
dsc_plots$response_fac <- as.factor(dsc_plots$response)

saveRDS(dsc_plots, paste0("../../output/sims_paper_figures_inter/", input, "_for_plotting.rds"))

