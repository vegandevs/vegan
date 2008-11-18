"print.vectorfit" <-
    function (x, ...) 
{
    if (x$permutations) 
        eps <- 1/x$permutations
    else eps <- .Machine$double.eps
    out <- cbind(x$arrows, r2 = x$r, "Pr(>r)" = x$pvals)
    printCoefmat(out, eps.Pvalue = eps, na.print = "", ...)
    if (x$permutations) {
        cat("P values based on", x$permutations, "permutations")
        if (!is.null(x$strata)) 
            cat(", stratified within", x$strata)
        cat(".\n")
    }
    invisible(x)
}
