`metaMDSiter` <-
    function (dist, k = 2, try = 20, trymax = 20, trace = 1, plot = FALSE,
              previous.best, engine = "monoMDS", maxit = 200,
              parallel = getOption("mc.cores"), ...)
{
    engine <- match.arg(engine, c("monoMDS", "isoMDS"))
    EPS <- 0.05
    if (engine == "monoMDS")
        EPS <- EPS/100 # monoMDS stress (0,1), isoMDS (0,100)
    RESLIM <- 0.01
    RMSELIM <- 0.005
    converged <- 0
    ## set tracing for engines
    isotrace <- max(0, trace - 1)
    monotrace <- engine == "monoMDS" && trace > 1
    ## explain monoMDS convergence codes (sol$icause)
    monomsg <- c("no. of iterations >= maxit",
                 "stress < smin",
                 "stress ratio > sratmax",
                 "scale factor of the gradient < sfgrmin")
    ## monoMDS trace >= 2
    monostop <- function(mod) {
        if (mod$maxits == 0)
            return(NULL)
        lab <- monomsg[mod$icause]
        cat("   ", mod$iters, "iterations: ", lab, "\n")
    }
    ## collect monoMDS convergence code for trace
    if (trace && engine == "monoMDS")
        stopcoz <- numeric(4)
    ## Previous best or initial configuration
    if (!missing(previous.best) && !is.null(previous.best)) {
        ## check if previous.best is from metaMDS or isoMDS
        if (inherits(previous.best, c("metaMDS", "monoMDS")) ||
            is.list(previous.best) &&
            all(c("points", "stress") %in% names(previous.best))) {
            ## Previous best may come from another 'engine' or
            ## 'model': extract its 'points' and use as an initial
            ## configuration with 'maxit = 0' to evaluate the stress
            ## in current case, or take a matrix as configuration.
            init <- previous.best$points
            bestry <- max(0, previous.best$bestry)
            trybase <- max(0, previous.best$tries)
            converged <- max(0, previous.best$converged)
            nc <- NCOL(init)
            if (nc > k)
                init <- init[, 1:k, drop = FALSE]
            else if (nc < k)
                for (i in 1:(k-nc))
                    init <- cbind(init, runif(NROW(init), -0.1, 0.1))
            if (trace)
                cat(sprintf("Starting from %d-dimensional configuration\n",
                            nc))
        } else {
            init <- as.matrix(previous.best)
            bestry <- 0
            trybase <- 0
        }
        ## evaluate stress
        s0 <- switch(engine,
                     "monoMDS" = monoMDS(dist, y = init, k = k, maxit = 0, ...),
                     "isoMDS" = isoMDS(dist, y = init, k = k, maxit = 0))
        ## Check whether model changed
        if (is.list(previous.best) && !is.null(previous.best$stress) &&
            !isTRUE(all.equal(previous.best$stress, s0$stress))) {
            if (trace) cat("Stress differs from 'previous.best': reset tries\n")
            if (inherits(previous.best, "metaMDS"))
                previous.best$tries <- 0
        }
    } else {
        ## no previous.best: start with cmdscale
        s0 <- switch(engine,
                 "monoMDS" = monoMDS(dist, y = cmdscale(dist, k = k), k = k,
                 maxit = maxit, ...),
                 "isoMDS" = isoMDS(dist, k = k, trace = isotrace,
                                   maxit = maxit))
        bestry <- 0
        trybase <- 0
    }
    if (trace)
        cat("Run 0 stress", s0$stress, "\n")
    if (monotrace)
        monostop(s0)
    tries <- 0
    ## Prepare for parallel processing
    if (is.null(parallel))
        parallel <- 1
    hasClus <- inherits(parallel, "cluster")
    isParal <- hasClus || parallel > 1
    isMulticore <- .Platform$OS.type == "unix" && !hasClus
    if (isParal && !isMulticore && !hasClus) {
        parallel <- makeCluster(parallel)
        clusterEvalQ(parallel, library(vegan))
    }
    ## get the number of clusters
    if (inherits(parallel, "cluster"))
        nclus <- length(parallel)
    else
        nclus <- parallel
    ## proper iterations
    while(tries < try || tries < trymax && converged == 0) {
        init <- replicate(nclus, initMDS(dist, k = k))
        if (nclus > 1) isotrace <- FALSE
        if (isParal) {
            if (isMulticore) {
                stry <-
                    mclapply(1:nclus, function(i)
                             switch(engine,
                                    "monoMDS" = monoMDS(dist, init[,,i], k = k,
                                    maxit = maxit, ...),
                                    "isoMDS" = isoMDS(dist, init[,,i], k = k,
                                    maxit = maxit, tol = 1e-07,
                                    trace = isotrace)),
                             mc.cores = parallel)
            } else {
                stry <-
                    parLapply(parallel, 1:nclus, function(i)
                              switch(engine,
                                     "monoMDS" = monoMDS(dist, init[,,i], k = k,
                                     maxit = maxit, ...),
                                     "isoMDS" = isoMDS(dist, init[,,i], k = k,
                                     maxit = maxit, tol = 1e-07, trace = isotrace)))
            }
        } else {
            stry <- list(switch(engine,
                                "monoMDS" = monoMDS(dist, init[,,1], k = k,
                                maxit = maxit, ...),
                                "isoMDS" = isoMDS(dist, init[,,1], k = k,
                                maxit = maxit, tol = 1e-07, trace = isotrace)))
        }
        ## analyse results of 'nclus' tries
        for (i in 1:nclus) {
            tries <- tries + 1
            if (trace)
                cat("Run", tries, "stress", stry[[i]]$stress, "\n")
            if (trace && engine == "monoMDS")
                stopcoz[stry[[i]]$icause] <- stopcoz[stry[[i]]$icause] + 1L
            if (monotrace)
                monostop(stry[[i]])
            if ((s0$stress - stry[[i]]$stress) > -EPS) {
                pro <- procrustes(s0, stry[[i]], symmetric = TRUE)
                if (plot && k > 1)
                    plot(pro)
                if (stry[[i]]$stress < s0$stress) {
                    s0 <- stry[[i]]
                    ## New best solution has not converged unless
                    ## proved later
                    converged <- 0
                    bestry <- tries + trybase
                    if (trace)
                        cat("... New best solution\n")
                }
                summ <- summary(pro)
                if (trace)
                    cat("... Procrustes: rmse", summ$rmse, " max resid",
                        max(summ$resid), "\n")
                if (summ$rmse < RMSELIM && max(summ$resid) < RESLIM) {
                    if (trace)
                        cat("... Similar to previous best\n")
                    converged <- converged + 1
                }
            }
            flush.console()
        }
    }
    if (trace) {
        if (converged > 0)
            cat("*** Best solution repeated", converged, "times\n")
        else if (engine == "monoMDS") {
            cat(sprintf(
                "*** Best solution was not repeated -- %s stopping criteria:\n",
                engine))
            for (i in seq_along(stopcoz))
                if (stopcoz[i] > 0)
                    cat(sprintf("%6d: %s\n", stopcoz[i], monomsg[i]))
        }
    }
    ## stop socket cluster
    if (isParal && !isMulticore && !hasClus)
        stopCluster(parallel)
    if (!missing(previous.best) && inherits(previous.best, "metaMDS")) {
        tries <- tries + previous.best$tries
    }
    out <- s0
    out$ndim = k
    out$data <- attr(dist, "commname")
    out$distance <- attr(dist, "method")
    out$converged <- converged
    out$tries <- tries
    out$bestry <- bestry
    out$engine <- engine
    out
}
