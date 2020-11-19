options(stringsAsFactors = FALSE)

###Load libraries
library(foreach)
library(doMC)
library(optparse)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--n_cores"), type="integer")
outparse <- parse_args(parser)


###Register cores
registerDoMC(outparse$n_cores)

###Load complete manifest
manifest <- read.table("../data/gtex-v8-manifest-full-X.txt", header=FALSE)

###Tissues to keep
to_keep <- c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Artery_Tibial", "Lung", 
             "Muscle_Skeletal", "Nerve_Tibial", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg", 
             "Thyroid", "Whole_Blood")

###Obtain missing Y rate
res <- foreach(i = 1:nrow(manifest), .combine='c') %dopar% {
  Y <- readRDS(manifest[i, 1])$y_res
  Y <- Y[, which(colnames(Y) %in% to_keep)]
  missing_rate <- sum(is.na(Y))/length(Y)
  missing_rate
}

###Filter manifest
manifest_missing <- cbind(manifest, res)

###Write out results
write.table(manifest_missing, "../data/gtex-v8-manifest-full-X-10tissues-missing-rate.txt", sep="\t", row.names = FALSE, col.names=FALSE, quote=FALSE)

