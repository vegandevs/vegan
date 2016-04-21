`as.mcmc.list.oecosimu` <-
    function(x)
{
    x <- x$oecosimu$simulated
    chains <- attr(x, "chains")
    if (is.null(chains))
        chains <- 1
    nsim <- dim(x)[2]
    niter <- nsim / chains
    x <- lapply(1:chains, function(i) {
        z <- x[,((i-1) * niter + 1):(i * niter)]
        attr(z, "mcpar") <- c(attr(x, "burnin") + attr(x, "thin"),
            attr(x, "burnin") + attr(x, "thin") * niter, attr(x, "thin"))
        attr(z, "class") <- "mcmc"
        z
    })
    class(x) <- c("mcmc.list", class(x))
    x
}
