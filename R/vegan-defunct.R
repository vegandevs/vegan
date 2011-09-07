## "new" permutation code was moved to package 'permute' in R 2.0-0.
## Here we list as defunct those functions that are not in 'permute'.

`permuted.index2` <- function (n, control = permControl()) 
    .Defunct("permute::shuffle", package="vegan")

`getNumObs` <- function(object, ...) 
    .Defunct("nobs", package = "vegan")
