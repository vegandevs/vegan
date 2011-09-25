`oecosimu` <-
    function(comm, nestfun, method, nsimul=99,
             burnin=0, thin=1, statistic = "statistic",
             alternative = c("two.sided", "less", "greater"),
             ...)
{
    alternative <- match.arg(alternative)
    nestfun <- match.fun(nestfun)
    applynestfun <-
        function(x, fun = nestfun, statistic = "statistic", ...) {
            tmp <- fun(x, ...)
            if (is.list(tmp))
                tmp[[statistic]]
            else
                tmp
    }

    nm <- nullmodel(comm, method)
    if (nm$commsim$binary)
        comm <- ifelse(comm > 0, 1L, 0L)
    
    ind <- nestfun(comm, ...)
    indstat <-
        if (is.list(ind))
            ind[[statistic]]
        else
            ind
    
    ## estimate thinning for "tswap" (trial swap)
    if (nm$commsim$method == "tswap") {
        checkbrd <-sum(designdist(comm, "(J-A)*(J-B)", 
                                  "binary"))
        M <- nm$ncol
        N <- nm$nrow
        checkbrd <- M * (M - 1) * N * (N - 1)/4/checkbrd
        thin <- round(thin * checkbrd)
        burnin <- round(burnin * checkbrd)
    }
    x <- simulate(nm, nsim = nsimul, burnin = burnin, thin = thin)

    simind <- apply(x, 3, applynestfun, fun = nestfun, statistic = statistic, ...) 
    simind <- matrix(simind, ncol = nsimul)

    if (nm$commsim$isSeq) {
        if (thin > 1)
            attr(simind, "thin") <- thin
        if (burnin > 0)
            attr(simind, "burnin") <- burnin
    }
    
    sd <- apply(simind, 1, sd, na.rm = TRUE)
    z <- (indstat - rowMeans(simind, na.rm = TRUE))/sd
    if (any(sd < sqrt(.Machine$double.eps)))
        z[sd < sqrt(.Machine$double.eps)] <- 0
    pless <- rowSums(indstat <= simind, na.rm = TRUE)
    pmore <- rowSums(indstat >= simind, na.rm = TRUE)
    if (any(is.na(simind))) {
        warning("some simulated values were NA and were removed")
        nsimul <- nsimul - rowSums(is.na(simind))
    }
    p <- switch(alternative,
                two.sided = 2*pmin(pless, pmore),
                less = pless,
                greater = pmore)
    p <- pmin(1, (p + 1)/(nsimul + 1))
    
    ## ADDITION: if z is NA then it is not correct to calculate p values
    ## try e.g. oecosimu(dune, sum, "permat")
    if (any(is.na(z)))
        p[is.na(z)] <- NA

    if (is.null(names(indstat)))
        names(indstat) <- statistic
    if (!is.list(ind))
        ind <- list(statistic = ind)
    ind$oecosimu <- list(z = z, pval = p, simulated=simind, method=nm$commsim$method,
                         statistic = indstat, alternative = alternative)
    class(ind) <- c("oecosimu", class(ind))
    ind
}

