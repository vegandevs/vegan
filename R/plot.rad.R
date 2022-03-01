"plot.rad" <-
    function(x, xlab="Rank", ylab="Abundance", log = "y", ...)
{
    rnk <- seq(along=x)
    plot(rnk, x, log=log, xlab=xlab, ylab=ylab, ...)
    out <- list(species = cbind(rnk, x))
    class(out) <- "ordiplot"
    invisible(out)
}
