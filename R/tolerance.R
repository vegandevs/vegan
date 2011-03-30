##' S3 generic for function to compute tolerances
##'
##' Brought this in here from analogue because of tolerance.cca
##'
##' @param x an R object
##' @param ... arguments passed to other methods
`tolerance` <- function(x, ...)
    UseMethod("tolerance")
