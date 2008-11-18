"plot.radline" <-
    function (x, xlab = "Rank", ylab = "Abundance", type = "b", ...) 
{
    rad <- x$y
    fit <- fitted(x)
    rnk <- seq(along = rad)
    plot(rnk, rad, log = "y", xlab = xlab, ylab = ylab, type = "n", 
         ...)
    out <- NULL
    if (type == "b" || type == "p") 
        out <- points(x, ...)
    if (type == "b" || type == "l") 
        lines(x, ...)
    invisible(out)
}
