#### Code from Gao Wang.
#### Prepare summary statistics 
#### for ED step.

###Set options
options(stringsAsFactors = FALSE)

###Set seed
set.seed(123)

###Load libraries
library(optparse)

###Parse arguments
parser <- OptionParser()
parser <- add_option(parser, c("--input_path"), type="character")
parser <- add_option(parser, c("--output"), type="character")
parser <- add_option(parser, c("--n_random"), type="integer", default=4)
outparse <- parse_args(parser)

input_path <- outparse$input_path
output <- outparse$output
num_random <- outparse$n_random

###Functions required for the actual step
matxMax <- function(mtx) {
  max_idx <- which.max(mtx)
  colmn <- max_idx %/% nrow(mtx) + 1
  row <- max_idx %% nrow(mtx)
  return( matrix(c(row, colmn), 1))
}

remove_rownames = function(x) {
  for (name in names(x)) rownames(x[[name]]) = NULL
  return(x)
}

extract_one_data = function(infile, n_random) {
  # If cannot read the input for some reason then let it go. I dont care losing one.
  dat = tryCatch(readRDS(infile)$sumstats, error = function(e) return(NULL))
  if (is.null(dat)) return(NULL)
  z = abs(dat$bhat/dat$shat)
  max_idx = matxMax(z)
  strong = list(bhat = dat$bhat[max_idx[1],,drop=FALSE], shat= dat$shat[max_idx[1],,drop=FALSE])
  if (max_idx[1] == 1) {
    sample_idx = 2:nrow(z)
  } else if (max_idx[1] == nrow(z)) {
    sample_idx = 1:(max_idx[1]-1)
  } else {
    sample_idx = c(1:(max_idx[1]-1), (max_idx[1]+1):nrow(z))
  }
  random_idx = sample(sample_idx, n_random, replace = TRUE)
  random = list(bhat = dat$bhat[random_idx,,drop=FALSE], shat = dat$shat[random_idx,,drop=FALSE])
  return(list(random = remove_rownames(random),  strong = remove_rownames(strong)))
}

merge_data = function(res, one_data) {
  if (length(res) == 0) {
    return(one_data)
  } else if (is.null(one_data)) {
    return(res)
  } else {
    for (d in names(one_data)) {
      for (s in names(one_data[[d]])) {
        res[[d]][[s]] = rbind(res[[d]][[s]], one_data[[d]][[s]])
      }
    }
    return(res)
  }
}


###Build strong and random effects
res = list()
for (f in list.files(path=input_path, pattern = "\\.rds$")) {
  res = merge_data(res, extract_one_data(paste0(input_path, "/", f), num_random))
}

###Reshape data to be in the format expected by the GTEx pipeline
res_pipe <- list(strong.b=res$strong$bhat, strong.s=res$strong$shat,
                 random.b=res$random$bhat, random.s=res$random$shat)

###Save results to a file
saveRDS(res_pipe, output)