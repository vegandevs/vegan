"print.summary.bioenv" <-
    function(x, ...)
{
    out <- data.frame(size = x$size, correlation = x$cor)
    rownames(out) <- x$var
    printCoefmat(out, ...)
    invisible(x)
}
