`oecosimu` <-
    function(comm, nestfun, method, nsimul=99,
             burnin=0, thin=1, statistic = "statistic",
             alternative = c("two.sided", "less", "greater"),
             parallel = getOption("mc.cores"), ...)
{
    alternative <- match.arg(alternative)
    nestfun <- match.fun(nestfun)
    if (length(statistic) > 1)
        stop("only one 'statistic' is allowed")
    applynestfun <-
        function(x, fun = nestfun, statistic = "statistic", ...) {
            tmp <- fun(x, ...)
            if (is.list(tmp))
                tmp[[statistic]]
            else
                tmp
    }
    if (inherits(comm, "simmat")) {
        x <- comm
        method <- attr(x, "method")
        nsimul <- dim(x)[3]
        if (nsimul == 1)
            stop("only one simulation in ", sQuote(deparse(substitute(comm))))
        comm <- attr(comm, "data")
        simmat_in <- TRUE
    } else {
        simmat_in <- FALSE
        if (inherits(comm, "nullmodel")) {
            nm <- comm
            comm <- comm$data
        } else {
            nm <- nullmodel(comm, method)
            if (nm$commsim$binary)
                comm <- nm$data
        }
        method <- nm$commsim$method
    }
    
    ind <- nestfun(comm, ...)
    indstat <-
        if (is.list(ind))
            ind[[statistic]]
        else
            ind
    if (!simmat_in) {
        ## estimate thinning for "tswap" (trial swap)
        if (method == "tswap") {
            checkbrd <-sum(designdist(comm, "(J-A)*(J-B)", 
                                      "binary"))
            M <- nm$ncol
            N <- nm$nrow
            checkbrd <- M * (M - 1) * N * (N - 1)/4/checkbrd
            thin <- round(thin * checkbrd)
            burnin <- round(burnin * checkbrd)
        }
        x <- simulate(nm, nsim = nsimul, burnin = burnin, thin = thin)
    }

    ## Go to parallel processing if 'parallel > 1' or 'parallel' could
    ## be a pre-defined socket cluster or 'parallel = NULL' in which
    ## case it could be setDefaultCluster (or a user error)
    if (is.null(parallel) && getRversion() >= "2.15.0")
        parallel <- get("default", envir = parallel:::.reg)
    ## still NULL? was not a defaultCluster, so make to non-parallel
    if (is.null(parallel) || getRversion() < "2.14.0")
        parallel <- 1
    hasClus <- inherits(parallel, "cluster")
    if ((hasClus || parallel > 1)  && require(parallel)) {
        if(.Platform$OS.type == "unix" && !hasClus) {
            tmp <- mclapply(1:nsimul,
                            function(i)
                            applynestfun(x[,,i], fun=nestfun,
                                         statistic = statistic, ...),
                            mc.cores = parallel)
            simind <- do.call(cbind, tmp)
        } else {
            ## if hasClus, do not set up and stop a temporary cluster
            if (!hasClus) {
                parallel <- makeCluster(parallel)
                ## make vegan functions available: others may be unavailable
                clusterEvalQ(parallel, library(vegan))
            }
            simind <- parApply(parallel, x, 3, function(z)
                               applynestfun(z, fun = nestfun,
                                            statistic = statistic, ...))
            if (!hasClus)
                stopCluster(parallel)
        }
    } else {
        simind <- apply(x, 3, applynestfun, fun = nestfun,
                        statistic = statistic, ...)
    }
    
    simind <- matrix(simind, ncol = nsimul)

    if (attr(x, "isSeq")) {
        attr(simind, "thin") <- attr(x, "thin")
        attr(simind, "burnin") <- attr(x, "start") - 1L
    }
    
    sd <- apply(simind, 1, sd, na.rm = TRUE)
    means <- rowMeans(simind, na.rm = TRUE)
    z <- (indstat - means)/sd
    if (any(sd < sqrt(.Machine$double.eps)))
        z[sd < sqrt(.Machine$double.eps)] <- 0
    pless <- rowSums(indstat >= simind, na.rm = TRUE)
    pmore <- rowSums(indstat <= simind, na.rm = TRUE)
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

    if (is.null(names(indstat)) && length(indstat) == 1)
        names(indstat) <- statistic
    oecosimu <- list(z = z, means = means, pval = p, simulated=simind,
                     method=method, statistic = indstat,
                     alternative = alternative, isSeq = attr(x, "isSeq"))
    out <- list(statistic = ind, oecosimu = oecosimu)
    attr(out, "call") <- match.call()
    class(out) <- "oecosimu"
    out
}

