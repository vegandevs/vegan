"plot.rad" <-
    function(x, xlab="Rank", ylab="Abundance",  ...)
{
    rnk <- seq(along=x)
    plot(rnk, x, log="y", xlab=xlab, ylab=ylab, ...)
    out <- list(species = cbind(rnk, x))
    class(out) <- "ordiplot"
    invisible(out)
}
