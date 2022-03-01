print.commsim <- function(x, ...) {
    cat("An object of class", dQuote(class(x)[1L]), "\n")
    isSeq <- ifelse(x$isSeq, "sequential", "non-sequential")
    if(x$binary)
        kind <- "binary"
    else
        kind <- ifelse(x$mode == "integer", "count", "abundance")
    cat(sQuote(x$method), " method (", 
        kind, ", ", isSeq, ", ", x$mode, " mode)\n\n", sep="")
    invisible(x)
}
