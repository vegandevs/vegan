`metaMDSiter` <-
    function (dist, k = 2, trymax = 20, trace = 1, plot = FALSE, 
              previous.best, engine = "monoMDS", ...) 
{
    engine <- match.arg(engine, c("monoMDS", "isoMDS"))
    if (engine == "isoMDS")
        require(MASS) || stop("Needs package MASS (function isoMDS)")
    EPS <- 0.05
    if (engine == "monoMDS")
        EPS <- EPS/100 # monoMDS stress (0,1), isoMDS (0,100) 
    RESLIM <- 0.01
    RMSELIM <- 0.005
    SOL <- FALSE
    converged <- FALSE
    isotrace <- max(0, trace - 1)
    ## Previous best or initial configuration 
    if (!missing(previous.best) && !is.null(previous.best)) {
        ## check if previous.best is from metaMDS or isoMDS
        if (inherits(previous.best, "metaMDS") ||
            is.list(previous.best) &&
            all(c("points", "stress") %in% names(previous.best))) {
            ## Previous best may come from another 'engine' or
            ## 'model': extract its 'points' and use as an initial
            ## configuration with 'maxit = 0' to evaluate the stress
            ## in current case, or take a matrix as configuration.
            init <- previous.best$points
            nc <- NCOL(init)
            if (nc > k)
                init <- init[, 1:k, drop = FALSE]
            else if (nc < k)
                for (i in 1:(k-nc))
                    init <- cbind(init, runif(NROW(init), -0.1, 0.1))
            if (trace)
                cat(gettextf("Starting from %d-dimensional configuration\n", nc))
        } else {
            init <- as.matrix(previous.best)
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
                 "monoMDS" = monoMDS(dist, y = cmdscale(dist, k = k), k = k, ...),
                 "isoMDS" = isoMDS(dist, k = k, trace = isotrace))
    }
    if (trace) 
        cat("Run 0 stress", s0$stress, "\n")
    tries <- 0
    while(tries < trymax) {
        tries <- tries + 1
        stry <- switch(engine,
                       "monoMDS" = monoMDS(dist, k = k, maxit = 200, ...),
                       "isoMDS" = isoMDS(dist, initMDS(dist, k = k), k = k,
                       maxit = 200, tol = 1e-07, trace = isotrace))
        if (trace) {
            cat("Run", tries, "stress", stry$stress, "\n")
        }
        if ((s0$stress - stry$stress) > -EPS) {
            pro <- procrustes(s0, stry, symmetric = TRUE)
            if (plot && k > 1) 
                plot(pro)
            if (stry$stress < s0$stress) {
                s0 <- stry
                if (trace) 
                    cat("... New best solution\n")
            }
            summ <- summary(pro)
            if (trace) 
                cat("... procrustes: rmse", summ$rmse, " max resid", 
                    max(summ$resid), "\n")
            if (summ$rmse < RMSELIM && max(summ$resid) < RESLIM) {
                if (trace) 
                    cat("*** Solution reached\n\n")
                converged <- TRUE
                break
            }
        }
        flush.console()
    }
    if (!missing(previous.best) && inherits(previous.best, "metaMDS")) {
        tries <- tries + previous.best$tries
    }
    out <- s0
    out$ndim = k
    out$data <- attr(dist, "commname")
    out$distance <- attr(dist, "method")
    out$converged <- converged
    out$tries <- tries
    out$engine <- engine
    out
}
