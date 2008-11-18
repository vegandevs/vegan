`print.adonis` <-
    function(x, ...)
{
    cat("\nCall:\n")
    cat(deparse(x$call), "\n\n")
    printCoefmat(x$aov.tab, eps = 1/nrow(x$f.perms), na.print = "")
    invisible(x)
}
