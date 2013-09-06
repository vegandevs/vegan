`print.vectorfit` <-
    function (x, ...) 
{
    out <- cbind(x$arrows, r2 = x$r, "Pr(>r)" = x$pvals)
    printCoefmat(out, na.print = "",
                 zap.ind = seq_len(ncol(out)-2), ...)
    if (x$permutations) {
        cat("P values based on", x$permutations, "permutations")
        if (!is.null(x$strata)) 
            cat(", stratified within", x$strata)
        cat(".\n")
    }
    invisible(x)
}
