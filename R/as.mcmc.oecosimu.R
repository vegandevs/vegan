`as.mcmc.oecosimu` <-
    function(x)
{
    chains <- attr(x$oecosimu$simulated, "chains")
    if (!is.null(chains) && chains > 1)
        stop("try as.mcmc.list for multiple chain")
    x <- as.ts(x)
    mcpar <- attr(x, "tsp")
    mcpar[3] <- round(1/mcpar[3])
    attr(x, "mcpar") <- mcpar
    class(x) <- c("mcmc", class(x))
    x
}
