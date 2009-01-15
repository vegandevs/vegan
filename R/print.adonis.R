`print.adonis` <-
    function(x, ...)
{
    cat("\nCall:\n")
    cat(deparse(x$call), "\n\n")
    printCoefmat(x$aov.tab,  na.print = "")
    invisible(x)
}
