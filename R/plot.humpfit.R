"plot.humpfit" <-
    function(x, xlab="Biomass", ylab="Species Richness", lwd=2, l.col="blue",
             p.col = 1, type="b", ...)
{
    plot(x$x, x$y, xlab = xlab, ylab = ylab, type="n", ...)
    if (type == "b" || type == "p")
        points(x, col = p.col, ...)
    if (type == "b" || type == "l")
        lines(x, lwd = lwd, col = l.col, ...)
    invisible()
}
