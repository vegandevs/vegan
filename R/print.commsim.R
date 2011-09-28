print.commsim <- function(x, ...) {
    cat("An object of class", dQuote(class(x)[1L]), "\n")
    isSeq <- ifelse(x$isSeq, "sequential", "non-sequential")
    binary <- ifelse(x$binary, "binary", "count")
    cat(sQuote(x$method), " method (", 
        binary, ", ", isSeq, ", ", x$mode, " mode)\n\n", sep="")
    invisible(x)
}
