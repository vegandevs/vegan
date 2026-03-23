`plot.rad` <-
    function(x, xlab="Rank", ylab="Abundance", log = "y", ...)
{
    rnk <- seq(along=x)
    plot(rnk, x, log=log, xlab=xlab, ylab=ylab, ...)
    out <- list(species = cbind(rnk, x))
    class(out) <- "ordiplot"
    invisible(out)
}

### plot.rad.frame could be implemented as a loop over rad models

`plot.rad.frame` <-
    function(x, ...)
{
    .NotYetImplemented()
}
