#### Sample individuals 
#### to be test set

###Set options
options(stringsAsFactors = FALSE)

###Set seed
set.seed(123)

###Load libraries
library(optparse)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--prop_testset"), type="numeric")
parser <- add_option(parser, c("--n"), type="integer")
parser <- add_option(parser, c("--output"), type="character")
outparse <- parse_args(parser)

prop_testset <- outparse$prop_testset
n <- outparse$n
output <- outparse$output

###Sample indeces for test set individuals
test_set <- sort(sample(x=c(1:n), size=round(n*prop_testset), replace=FALSE))

###Save index to file
saveRDS(test_set, output)