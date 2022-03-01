`points.radline` <-
    function (x, ...) 
{
    poi <- x$y
    rnk <- seq(along = poi)
    points(rnk, poi, ...)
    out <- list(species = cbind(rnk, poi))
    class(out) <- "ordiplot"
    invisible(out)
}

`points.radfit` <-
    function(x, ...)
{
    points.radline(x, ...)
}
