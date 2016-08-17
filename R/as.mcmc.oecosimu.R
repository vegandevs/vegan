`as.mcmc.oecosimu` <-
    function(x)
{
    ## mcmc only for sequential methods
    if (!x$oecosimu$isSeq)
        stop("as.mcmc available only for sequential null models")
    ## named variables
    rownames(x$oecosimu$simulated) <- names(x$oecosimu$z)
    chains <- attr(x$oecosimu$simulated, "chains")
    ## chains: will make each chain as an mcmc object and combine
    ## these to an mcmc.list
    if (!is.null(chains) && chains > 1) {
        x <- x$oecosimu$simulated
        nsim <- dim(x)[2]
        niter <- nsim / chains
        ## iterate over chains
        x <- lapply(1:chains, function(i) {
                        z <- x[, ((i-1) * niter + 1):(i * niter), drop = FALSE]
                        attr(z, "mcpar") <-
                            c(attr(x, "burnin") + attr(x, "thin"),
                              attr(x, "burnin") + attr(x, "thin") * niter,
                              attr(x, "thin"))
                        attr(z, "class") <- c("mcmc", class(z))
                        t(z)
                    })
        ## combine list of mcmc objects to a coda mcmc.list
        #x <- as.mcmc.list(x)
        class(x) <- "mcmc.list"
    } else { # one chain: make to a single mcmc object
        x <- as.ts(x)
        mcpar <- attr(x, "tsp")
        mcpar[3] <- round(1/mcpar[3])
        attr(x, "mcpar") <- mcpar
        class(x) <- c("mcmc", class(x))
    }
    x
}
