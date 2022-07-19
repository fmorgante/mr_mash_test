options(stringsAsFactors = FALSE)

dat <- read.table("../data/gtex-v8-full-X-number-of-snps.txt", header=TRUE, sep="\t")

res <- do.call(cbind, lapply(dat[,-1], summary))

write.table(res, "../data/gtex-v8-full-X-number-of-snps-summary-stats.txt", sep="\t", row.names = TRUE, col.names=TRUE, quote=FALSE)
