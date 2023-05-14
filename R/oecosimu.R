`oecosimu` <-
    function(comm, nestfun, method, nsimul=99,
             burnin=0, thin=1, statistic = "statistic",
             alternative = c("two.sided", "less", "greater"),
             batchsize = NA,
             parallel = getOption("mc.cores"), ...)
{
    alternative <- match.arg(alternative)
    nestfun <- match.fun(nestfun)
    if (length(statistic) > 1)
        stop("only one 'statistic' is allowed")
    if (!is.na(batchsize))
        batchsize <- batchsize * 1024 * 1024
    applynestfun <-
        function(x, fun = nestfun, statistic = "statistic", ...) {
            tmp <- fun(x, ...)
            if (is.list(tmp))
                tmp[[statistic]]
            else
                tmp
    }
    chains <- NULL
    if (inherits(comm, "simmat")) {
        x <- comm
        method <- attr(x, "method")
        nsimul <- dim(x)[3]
        if (nsimul == 1)
            stop(gettextf("only one simulation in '%s'",
                          deparse(substitute(comm))))
        comm <- attr(comm, "data")
        #thin <- attr(comm, "thin")
        burnin <- attr(x, "start") - attr(x, "thin")
        chains <- attr(x, "chains")
        simmat_in <- TRUE
    } else {
        simmat_in <- FALSE
        if (inherits(comm, "nullmodel")) {
            nm <- comm
            comm <- comm$data
        } else {
            nm <- nullmodel(comm, method)
            if (nm$commsim$binary) {
                ## sometimes people do not realize that null model
                ## makes their data binary
                if (max(abs(comm - nm$data)) > 0.1)
                    warning("nullmodel transformed 'comm' to binary data")
                comm <- nm$data
            }
        }
        method <- nm$commsim$method
    }
    ## Check the number of batches needed to run the requested number
    ## of simulations without exceeding arg 'batchsize', and find the
    ## size of each batch.
    if (!simmat_in && !is.na(batchsize)) {
        commsize <- object.size(comm)
        totsize <- commsize * nsimul
        if (totsize > batchsize) {
            nbatch <- ceiling(unclass(totsize/batchsize))
            batches <- diff(round(seq(0, nsimul, by = nsimul/nbatch)))
        } else {
            nbatch <- 1
        }
    } else {
        nbatch <- 1
    }
    if (nbatch == 1)
        batches <- nsimul

    ind <- nestfun(comm, ...)
    indstat <-
        if (is.list(ind))
            ind[[statistic]]
        else
            ind
    ## burnin of sequential models
    if (!simmat_in && nm$commsim$isSeq) {
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
        if (burnin > 0)
            nm <- update(nm, burnin)
    }
    ## start with empty simind
    simind <- NULL
    ## Go to parallel processing if 'parallel > 1' or 'parallel' could
    ## be a pre-defined socket cluster or 'parallel = NULL'.
    if (is.null(parallel))
        parallel <- 1
    hasClus <- inherits(parallel, "cluster")
    if (hasClus || parallel > 1) {
        if(.Platform$OS.type == "unix" && !hasClus) {
            for (i in seq_len(nbatch)) {
                ## simulate if no simmat_in
                if(!simmat_in)
                    x <- simulate(nm, nsim = batches[i], thin = thin)
                tmp <- mclapply(seq_len(batches[i]),
                                function(j)
                                applynestfun(x[,,j], fun=nestfun,
                                             statistic = statistic, ...),
                                mc.cores = parallel)
                simind <- cbind(simind, do.call(cbind, tmp))
            }
        } else {
            ## if hasClus, do not set up and stop a temporary cluster
            if (!hasClus) {
                parallel <- makeCluster(parallel)
                ## make vegan functions available: others may be unavailable
                clusterEvalQ(parallel, library(vegan))
            }
            for(i in seq_len(nbatch)) {
                if (!simmat_in)
                    x <- simulate(nm, nsim = batches[i], thin = thin)
                simind <- cbind(simind,
                                parApply(parallel, x, 3, function(z)
                                         applynestfun(z, fun = nestfun,
                                                      statistic = statistic, ...)))
            }
            if (!hasClus)
                stopCluster(parallel)
        }
    } else {
        for(i in seq_len(nbatch)) {
            ## do not simulate if x was already a simulation
            if(!simmat_in)
                x <- simulate(nm, nsim = batches[i], thin = thin)
            simind <- cbind(simind, apply(x, 3, applynestfun, fun = nestfun,
                                          statistic = statistic, ...))
        }
    }

    simind <- matrix(simind, ncol = nsimul)

    if (attr(x, "isSeq")) {
        attr(simind, "thin") <- attr(x, "thin")
        attr(simind, "burnin") <- burnin
        attr(simind, "chains") <- chains
    }

    sd <- apply(simind, 1, sd, na.rm = TRUE)
    means <- rowMeans(simind, na.rm = TRUE)
    z <- (indstat - means)/sd
    if (any(sd < sqrt(.Machine$double.eps)))
        z[sd < sqrt(.Machine$double.eps)] <- 0
    ## results can be integers or real: comparisons differ
    if (is.integer(indstat) && is.integer(simind)) {
        pless <- rowSums(indstat >= simind, na.rm = TRUE)
        pmore <- rowSums(indstat <= simind, na.rm = TRUE)
    } else {
        EPS <- sqrt(.Machine$double.eps)
        pless <- rowSums(indstat + EPS >= simind, na.rm = TRUE)
        pmore <- rowSums(indstat - EPS <= simind, na.rm = TRUE)
    }
    if (any(is.na(simind))) {
        warning("some simulated values were NA and were removed")
        nsimul <- nsimul - rowSums(is.na(simind))
    }
    p <- switch(alternative,
                two.sided = 2*pmin.int(pless, pmore),
                less = pless,
                greater = pmore)
    p <- pmin.int(1, (p + 1)/(nsimul + 1))

    ## ADDITION: if z is NA then it is not correct to calculate p values
    ## try e.g. oecosimu(dune, sum, "permat")
    if (any(is.na(z)))
        p[is.na(z)] <- NA

    ## take care that statistics have name, or some support functions
    ## can fail
    if (is.null(names(indstat))) {
        if (length(indstat) == 1)
            names(indstat) <- statistic
        else if (length(indstat) <= length(letters))
            names(indstat) <- letters[seq_along(indstat)]
        else
            names(indstat) <- paste0("stat", seq_along(indstat))
    }
    oecosimu <- list(z = z, means = means, pval = p, simulated=simind,
                     method=method, statistic = indstat,
                     alternative = alternative, isSeq = attr(x, "isSeq"))
    out <- list(statistic = ind, oecosimu = oecosimu)
    attr(out, "call") <- match.call()
    class(out) <- "oecosimu"
    out
}

