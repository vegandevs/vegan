`density.oecosimu` <-
    function(x, ...)
{
    out <- density(t(x$oecosimu$simulated), ...)
    out$call <- match.call()
    out$call[[1]] <- as.name("density")
    out
}

