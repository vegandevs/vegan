print.summary.clamtest <- function(x, digits=max(3, getOption("digits") - 3), ...) {
    cat("Two Groups Species Classification Method (CLAM)\n\n")
    cat("Specialization threshold =", x$specialization)
    cat("\nAlpha level =", x$alpha)
    cat("\n\nEstimated sample coverage:\n")
    print(x$coverage, digits=digits)
    cat("\nMinimum abundance for classification:\n")
    print(structure(c(x$minv[[1]][1,2], x$minv[[2]][1,1]),
        .Names=x$labels))
    cat("\n")
    printCoefmat(x$summary, digits=digits, ...)
}

