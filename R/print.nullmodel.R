print.nullmodel <- function(x, ...) {
    isSeq <- ifelse(x$commsim$isSeq, "sequential", "non-sequential")
    binary <- ifelse(x$commsim$binary, "binary", "count")
    cat("An object of class \"", class(x)[1L], "\"\n", sep="")
    cat("\"", x$commsim$method, "\" method (", 
        binary, ", ", isSeq, ")\n", sep="")
    cat(x$nrow, "x", x$ncol, "matrix\n")
    if (x$commsim$isSeq)
        cat("Iterations =", x$iter, "\n\n") else cat("\n")
    invisible(x)
}
