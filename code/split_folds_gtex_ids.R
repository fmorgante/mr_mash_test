###########################
### Fabio Morgante
### 05/12/2021
###
### Split GTEx individuals into
### 5 folds for cross-validation
###########################

options(stringsAsFactors = FALSE)
set.seed(1)

###Read in the data
gtex_ids <- read.table("../data/gtex-v8-ids.txt", header = FALSE, sep="\t")

###Split the data into 5-folds
gtex_ids_shuf <- gtex_ids[sample(nrow(gtex_ids)), , drop=FALSE]
folds <- cut(seq(1, nrow(gtex_ids_shuf)), breaks=5, labels=FALSE)
gtex_ids_shuf_folds <- cbind(gtex_ids_shuf, folds)
colnames(gtex_ids_shuf_folds) <- c("id", "fold")
gtex_ids_shuf_folds <- gtex_ids_shuf_folds[order(gtex_ids_shuf_folds$id),]

###Write out the new data
write.table(gtex_ids_shuf_folds, "../data/gtex-v8-ids-folds.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")


