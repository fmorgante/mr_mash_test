options(stringsAsFactors = FALSE)

###Load libraries
library(foreach)
library(doMC)
library(optparse)
library(dplyr)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--n_cores"), type="integer")
outparse <- parse_args(parser)


###Register cores
registerDoMC(outparse$n_cores)

###Load complete manifest
manifest <- read.table("../data/gtex-v8-manifest.txt", header=FALSE)

###Obtain gene names
gene_names <- vector("character", nrow(manifest))
for(j in 1:nrow(manifest)){
  gene_names[j] <- paste(head(unlist(strsplit(tail(unlist(strsplit(manifest[j, 1], "/")), 1), ".", fixed=TRUE)), 2), collapse=".")
}
manifest <- cbind(manifest, gene_names)

###Obtain number of non-missing samples
res <- foreach(i = 1:nrow(manifest), .combine='bind_rows') %dopar% {
  Y <- readRDS(manifest[i, 1])$y_res
  non_missing <- apply(Y, 2, function(x){sum(!is.na(x))})
  data.frame(Gene_Name=manifest[i, 2], rbind(non_missing))
}

###Write out results
write.table(gtex_non_missing, "../data/gtex-v8-non-missing-samples-by-tissue.txt", sep="\t", row.names = FALSE, col.names=TRUE, quote=FALSE)

