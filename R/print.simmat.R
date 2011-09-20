print.simmat <- function(x, ...) {
    isSeq <- ifelse(attr(x, "isSeq"), "sequential", "non-sequential")
    binary <- ifelse(attr(x, "binary"), "binary", "count")
    d <- dim(x)
    cat("An object of class \"", class(x)[1L], "\"\n", sep="")
    cat("\"", attr(x, "method"), "\" method (", 
        binary, ", ", isSeq, ")\n", sep="")
    cat(d[1L], "x", d[2L], "matrix\n")
    cat("Number of permuted matrices =", d[3L], "\n")
    if (attr(x, "isSeq")) {
        cat("Start = ", attr(x, "start"), ", End = ", attr(x, "end"), 
            ", Thin = ", attr(x, "thin"), "\n\n", sep="") 
        } else cat("\n")
    invisible(x)
}
