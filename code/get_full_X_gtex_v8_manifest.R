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
manifest <- read.table("../data/gtex-v8-manifest.txt", header=FALSE)

###Obtain indexes of genes with full X (i.e., 838 individuals)
res <- foreach(i = 1:nrow(manifest), .combine='c') %dopar% {
  X <- readRDS(manifest[i, 1])$X
  if(nrow(X)==838)
   i
}

###Filter manifest
manifest_filt <- manifest[res, ]

###Write out results
write.table(manifest_filt, "../data/gtex-v8-manifest-full-X.txt", sep="\t", row.names = FALSE, col.names=FALSE, quote=FALSE)
