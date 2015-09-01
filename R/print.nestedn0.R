"print.nestedn0" <-
function(x, ...)
{
    cat("Nestedness index N0:", format(x$statistic), "\n")
    invisible(x)
}

