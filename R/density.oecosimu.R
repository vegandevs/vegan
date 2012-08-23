`density.oecosimu` <-
    function(x, ...)
{
    cols <- nrow(x$oecosimu$simulated)
    if (cols > 1)
        warning("'density' is meaningful only with one statistic, you have ", cols)
    out <- density(t(x$oecosimu$simulated), ...)
    out$call <- match.call()
    out$call[[1]] <- as.name("density")
    out
}

