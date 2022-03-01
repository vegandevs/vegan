"print.nestedtemp" <-
function(x, ...)
{
    cat("nestedness temperature:", format(x$statistic, ...), "\n")
    cat("with matrix fill", format(x$fill, ...), "\n")
    invisible(x)
}

