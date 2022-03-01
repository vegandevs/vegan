print.nullmodel <- function(x, ...) {
    isSeq <- ifelse(x$commsim$isSeq, "sequential", "non-sequential")
    if (x$commsim$binary)
        kind <- "binary"
    else
        kind <- ifelse(x$commsim$mode == "integer", "count", "abundance")
    cat("An object of class", dQuote(class(x)[1L]), "\n")
    cat(sQuote(x$commsim$method), " method (", 
        kind, ", ", isSeq, ")\n", sep="")
    cat(x$nrow, "x", x$ncol, "matrix\n")
    if (x$commsim$isSeq)
        cat("Iterations =", x$iter, "\n\n") else cat("\n")
    invisible(x)
}
