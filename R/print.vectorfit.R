`print.vectorfit` <-
    function (x, ...) 
{
    out <- cbind(x$arrows, r2 = x$r, "Pr(>r)" = x$pvals)
    printCoefmat(out, na.print = "",
                 zap.ind = seq_len(ncol(out)-2), ...)
    if (x$permutations) {
        cat(howHead(x$control))
    }
    invisible(x)
}
