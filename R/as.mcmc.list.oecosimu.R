`as.mcmc.oecosimu` <-
    function(x)
{
    x <- x$oecosimu$simulated
    chains <- attr(x, "chains")
    if (is.null(chains))
        chains <- 1
    nsim <- dim(x)[2]
    niter <- nsim / chains
    x <- lapply(1:chains, function(i)
        as.mcmc(x[,((i-1) * niter + 1):(i * niter)]))
    class(x) <- c("mcmc.list", class(x))
    x
}


smfun <- function(x, burnin, nsim, thin) {
    nm <- nullmodel(x, "swap")
    nm <- update(nm, nsim=burnin)
    simulate(nm, nsim=nsim, thin=thin)
}
smlist <- replicate(3, smfun(x, burnin=50, nsim=10, thin=5), simplify=FALSE)
smbind(smlist, MARGIN=3) # Number of permuted matrices = 30

## need to check start/end/thin to be correct with summary.mcmc.list

