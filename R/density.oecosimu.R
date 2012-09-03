`density.oecosimu` <-
    function(x, ...)
{
    cols <- nrow(x$oecosimu$simulated)
    if (cols > 1)
        warning("'density' is meaningful only with one statistic, you have ", cols)
    obs <- x$oecosimu$statistic
    out <- density(rbind(obs, t(x$oecosimu$simulated)), ...)
    out$observed <- obs
    out$call <- match.call()
    out$call[[1]] <- as.name("density")
    class(out) <- c("vegandensity", class(out))
    out
}
