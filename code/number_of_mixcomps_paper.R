###Get number of data-driven matrices in GTEx application

fold1 <- readRDS("/Users/fmorgante/Desktop/fold_1.ted_unconstrained.rds")
fold2 <- readRDS("/Users/fmorgante/Desktop/fold_2.ted_unconstrained.rds")
fold3 <- readRDS("/Users/fmorgante/Desktop/fold_3.ted_unconstrained.rds")
fold4 <- readRDS("/Users/fmorgante/Desktop/fold_4.ted_unconstrained.rds")
fold5 <- readRDS("/Users/fmorgante/Desktop/fold_5.ted_unconstrained.rds")

c(length(fold1$U), length(fold2$U), length(fold3$U), length(fold4$U), length(fold5$U))


###Get number of data-driven matrices in simulations (GTExrealX_indepV_sharedB_3causalrespr10_indepReps)
data_driven <- vector("numeric", 20)

for(i in 1:20){
  dat <- readRDS(paste0("/project/mstephens/fmorgante/mr_mash_test/output/test/mvreg_all_genes_prior_GTExrealX_indepV_sharedB_3causalrespr10_indepReps_inter/rep", i, "/prior/matrices/mvreg_all_genes_prior_GTExrealX_indepV_sharedB_3causalrespr10_indepReps.EZ.prior.rds"))
  dat[c("identity", "Y1", "Y2", "Y3", "Y4", "Y5", "Y6", "Y7", "Y8", "Y9", "Y10", 
            "equal_effects", "simple_het_1", "simple_het_2", "simple_het_3")] <- NULL
  
  data_driven[i] <- length(dat)
}

summary(data_driven)

###Get number of data-driven matrices in simulations (GTExrealX_indepV_sharedB_r10_indepReps)
data_driven <- vector("numeric", 20)

for(i in 1:20){
  dat <- readRDS(paste0("/project/mstephens/fmorgante/mr_mash_test/output/test/mvreg_all_genes_prior_GTExrealX_indepV_sharedB_r10_indepReps_inter/rep", i, "/prior/matrices/mvreg_all_genes_prior_GTExrealX_indepV_sharedB_r10_indepReps.EZ.prior.rds"))
  dat[c("identity", "Y1", "Y2", "Y3", "Y4", "Y5", "Y6", "Y7", "Y8", "Y9", "Y10", 
        "equal_effects", "simple_het_1", "simple_het_2", "simple_het_3")] <- NULL
  
  data_driven[i] <- length(dat)
}

summary(data_driven)

