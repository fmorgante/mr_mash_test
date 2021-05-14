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

###Obtain size of Y
res <- foreach(i = 1:nrow(manifest), .combine='rbind') %dopar% {
  Y <- readRDS(manifest[i, 1])$y_res
  data.frame(gene=manifest[i, 1], ntissues=ncol(Y))
}

###Write out results
write.table(res, "../data/gtex-v8-manifest-tissue-number.txt", sep="\t", row.names = FALSE, col.names=FALSE, quote=FALSE)
