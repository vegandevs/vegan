"print.factorfit" <-
    function (x, ...) 
{
    cat("Centroids:\n")
    printCoefmat(x$centroids, tst.ind = 1:ncol(x$centroids), na.print = "", ...)
    cat("\nGoodness of fit:\n")
    out <- cbind(r2 = x$r, "Pr(>r)" = x$pvals)
    if (x$permutations) {
        printCoefmat(out, has.Pvalue = TRUE, ...)
        cat("P values based on", x$permutations, "permutations")
        if (!is.null(x$strata)) 
            cat(", stratified within", x$strata)
        cat(".\n")
    }
    else  printCoefmat(out, na.print = "", ...)
    invisible(x)
}
