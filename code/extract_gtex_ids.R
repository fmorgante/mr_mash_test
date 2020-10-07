options(stringsAsFactors = FALSE)

###Load filtered manifest
manifest <- read.table("../data/gtex-v8-manifest-full-X.txt", header=FALSE)

###Load first gene's X and extract the individual ID
X <- readRDS(manifest[1, 1])$X
ids <- rownames(X)

write.table(ids, "../data/gtex-v8-ids.txt", sep="\t", row.names = FALSE, col.names=FALSE, quote=FALSE)
