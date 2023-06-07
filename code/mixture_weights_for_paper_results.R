dat <- readRDS("/project/mstephens/fmorgante/mr_mash_test/output/test/mvreg_all_genes_prior_GTExrealX_indepV_sharedB_3causalrespr10_indepReps_inter/rep1/prior/matrices/mvreg_all_genes_prior_GTExrealX_indepV_sharedB_3causalrespr10_indepReps.EZ.prior.rds")
dat1 <- readRDS("/project/mstephens/fmorgante/mr_mash_test/output/test/mvreg_all_genes_prior_GTExrealX_indepV_sharedB_3causalrespr10_indepReps/rep1/mr_mash_em_data/extract_X_1_indepV_sharedB_subsetcausalr_1_mlasso_init_1_univ_sumstats_1_mr_mash_em_data_1.rds")


e <- rep(0,length(dat))
u <- c(rep(1,3),rep(0,7))
U <- tcrossprod(u)
names(e) <- names(dat)
for (i in names(dat)) {
  dat[[i]] <- dat[[i]]/max(dat[[i]])
  e[i] <- mean(abs(U - dat[[i]]))
}
print(data.frame(e = round(e,digits = 3)))


head(sort(dat1$fit_obj$w0,decreasing=TRUE), 20)

matrices_tokeep <- c("ED_XtX.", "ED_tFLASH_default.", "ED_FLASH_default_1.",
                     "ED_tFLASH_nonneg.", "ED_FLASH_nonneg_1.", "ED_PCA_1.", 
                     "ED_tPCA.")
tokeep <- stringr::str_detect(names(dat1$fit_obj$w0), stringr::str_c(matrices_tokeep, collapse = "|"))
sum(dat1$fit_obj$w0[tokeep])/(1 - dat1$fit_obj$w0["null"])


