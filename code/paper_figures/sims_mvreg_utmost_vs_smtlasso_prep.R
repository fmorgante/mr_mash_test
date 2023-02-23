###Set options
options(stringsAsFactors=FALSE)

###Load libraries
library(dscrutils)
library(optparse)

###Function to convert dscquery output from list to data.frame suitable for plotting
convert_dsc_to_dataframe <- function(dsc){
  ###Data.frame to store the results after convertion
  dsc_df <- data.frame()

  ###Get length of list elements
  n_elem <- length(dsc$DSC)

  ###Loop through the dsc list
  for(i in 1:n_elem){
    ##Prepare vectors making up the final data frame
    r_scalar <- dsc$simulate.r[i]
    repp <- rep(dsc$DSC[i], times=r_scalar)
    n <- rep(dsc$simulate.n[i], times=r_scalar)
    p <- rep(dsc$simulate.p[i], times=r_scalar)
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
    df <- data.frame(rep=repp, n=n, p=p, p_num_caus=p_causal, r=r, response=response, pve=pve,
                     scenario=simulate, method=fit, score_metric=score, score_value=score.err, time=timing)
    dsc_df <- rbind(dsc_df, df)
  }

  return(dsc_df)
}

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--input"), type="character")
outparse <- parse_args(parser)
input <- outparse$input


###Load data
dsc_out <- dscquery("../../output/mvreg_mrmash_vs_mtlasso_vs_utmost", 
		    		c("simulate.n", "simulate.p", "simulate.p_causal", "simulate.r",   
                      "simulate.pve", "simulate.B_cor", "simulate.B_scale",
                      "simulate.X_cor", "simulate.X_scale",
                      "simulate.V_cor", "simulate.prop_testset",
                      "simulate", "fit", "score", "score.err", "fit.time"), 
                    conditions = paste0("$(simulate) == '", input, "'"), verbose=FALSE)


###Convert dscquery output from list to data.frame suitable for plotting
dsc_plots <- convert_dsc_to_dataframe(dsc_out)


###Create factor version of method
dsc_plots$method_fac <- factor(dsc_plots$method, levels=c("mr_mash", "mtlasso", "utmost"),
                               labels=c("mr_mash_can", "mtlasso", "enet"))

###Create factor version of response
dsc_plots$response_fac <- as.factor(dsc_plots$response)

saveRDS(dsc_plots, paste0("../../output/sims_paper_figures_inter/mvreg_", input, "_utmost_smtlasso_for_plotting.rds"))
