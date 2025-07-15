`plot.varpart` <-
    function(x, Xnames = x$tables, ...)
{
    plot(x$part, Xnames = Xnames, ...)
    invisible()
}

