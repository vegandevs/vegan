`print.anosim` <-
    function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\nCall:\n")
    cat(deparse(x$call), "\n")
    cat("Dissimilarity:", x$dissimilarity,"\n\n")
    cat("ANOSIM statistic R: ")
    cat(formatC(x$statistic, digits = digits), "\n")
    nperm <- x$permutations
    if (nperm) {
        cat("      Significance:", format.pval(x$signif), 
            "\n\n")
        cat(howHead(x$control))
    }
    cat("\n")
    invisible(x)
}
