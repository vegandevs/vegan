`print.adonis` <-
    function(x, ...)
{
    cat("\nCall:\n")
    cat(deparse(x$call), "\n\n")
    print(x$aov.tab)
    invisible(x)
}
