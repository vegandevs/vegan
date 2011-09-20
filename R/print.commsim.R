print.commsim <- function(x, ...) {
    cat("An object of class \"", class(x)[1L] , "\"\n", sep="")
    isSeq <- ifelse(x$isSeq, "sequential", "non-sequential")
    binary <- ifelse(x$binary, "binary", "count")
    cat("\"", x$method, "\" method (", 
        binary, ", ", isSeq, ", ", x$mode, " mode)\n\n", sep="")
    invisible(x)
}
