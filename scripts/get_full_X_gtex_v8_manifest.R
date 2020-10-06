options(stringsAsFactors = FALSE)

###Load complete manifest
manifest <- read.table("../data/gtex-v8-manifest.txt", header=FALSE)

res <- c()

###Obtain indexes of genes with full X (i.e., 838 individuals)
for(i in 1:nrow(manifest)){
  X <- readRDS(manifest[i, 1])$X
  if(nrow(X)==838)
   res <- c(res, i)
  rm(X)
  cat("Finished gene", i, "\n")
}

###Filter manifest
manifest_filt <- manifest[res, ]

###Write out results
write.table(manifest_filt, "../data/gtex-v8-manifest-full-X.txt", sep="\t", row.names = FALSE, col.names=FALSE, quote=FALSE)
