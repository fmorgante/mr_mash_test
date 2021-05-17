###########################
### Fabio Morgante
### 05/17/2021
###
### Obtain list of genes expressed
### in more than 1 tissue
###########################

options(stringsAsFactors=FALSE)

###Load data
dat <- read.table("../data/gtex-v8-manifest-tissue-number.txt", header=FALSE, sep="\t") 

###Extract genes expressed in more than 1 tissue
tokeep <- which(dat[, 2] > 1)
dat_filt <- dat[tokeep, ]

###Write the results to a table
write.table(dat_filt, "../data/gtex-v8-manifest-2ormore-tissues.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
