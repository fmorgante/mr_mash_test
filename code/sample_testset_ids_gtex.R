#### Sample individuals 
#### to be test set

###Set options
options(stringsAsFactors = FALSE)

###Load libraries
library(optparse)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--n_testset"), type="integer")
parser <- add_option(parser, c("--input"), type="character")
parser <- add_option(parser, c("--output"), type="character")
parser <- add_option(parser, c("--seed"), type="integer")
outparse <- parse_args(parser)

n_testset <- outparse$n_testset
n <- outparse$n
input <- outparse$input
output <- outparse$output
ranseed <- outparse$seed

###Set seed
set.seed(ranseed)

###Load the data
ids <- read.table(input, header=FALSE, sep="\t")
ids <- ids[, 1]

###Sample indeces for test set individuals
test_set <- sort(sample(x=ids, size=n_testset, replace=FALSE))

###Save index to file
saveRDS(test_set, output)