"print.anosim" <-
    function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\nCall:\n")
    cat(deparse(x$call), "\n")
    cat("Dissimilarity:", x$dissimilarity,"\n\n")
    cat("ANOSIM statistic R: ")
    cat(formatC(x$statistic, digits = digits), "\n")
    nperm <- x$permutations
    if (nperm) {
        cat("      Significance:", format.pval(x$signif, eps = 1/nperm), 
            "\n\n")
        cat("Based on ", nperm, " permutations")
    }
    if (!is.null(x$strata)) 
        cat(", stratified within", x$strata)
    cat("\n\n")
    invisible(x)
}
