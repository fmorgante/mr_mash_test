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

###Obtain missing Y rate
res <- foreach(i = 1:nrow(manifest), .combine='bind_rows') %dopar% {
  Y <- readRDS(manifest[i, 1])$y_res
  missing_rate <- apply(Y, 2, function(x){sum(is.na(x))})/nrow(Y)
  data.frame(Gene_Name=manifest[i, 2], rbind(missing_rate))
}

###Write out results
write.table(gtex_missing_rate, "../data/gtex-v8-missing-rate-by-tissue.txt", sep="\t", row.names = FALSE, col.names=TRUE, quote=FALSE)

