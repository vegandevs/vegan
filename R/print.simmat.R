print.simmat <- function(x, ...) {
    isSeq <- ifelse(attr(x, "isSeq"), "sequential", "non-sequential")
    if (attr(x, "binary")) {
        kind <- "binary"
    } else {
        kind <- ifelse(attr(x, "mode") == "integer", "count", "abundance")
    }
    d <- dim(x)
    cat("An object of class", dQuote(class(x)[1L]), "\n")
    cat(sQuote(attr(x, "method")), " method (", 
        kind, ", ", isSeq, ")\n", sep="")
    cat(d[1L], "x", d[2L], "matrix\n")
    cat("Number of permuted matrices =", d[3L], "\n")
    if (attr(x, "isSeq")) {
        chainInfo <- ""
        if (!is.null(attr(x, "chains")) && attr(x, "chains") > 1L)
            chainInfo <- paste0(" (", attr(x, "chains"), " chains)")
        cat("Start = ", attr(x, "start"), ", End = ", attr(x, "end"), 
            ", Thin = ", attr(x, "thin"), chainInfo, "\n\n", sep="") 
        } else cat("\n")
    invisible(x)
}
