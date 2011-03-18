`metaMDSiter` <-
    function (dist, k = 2, trymax = 20, trace = 1, plot = FALSE, 
              previous.best, ...) 
{
    if (!require(MASS)) 
        stop("Needs package MASS (function isoMDS): not found")
    EPS <- 0.05
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
            ## real "previous best"
            if (NCOL(previous.best$points) == k) {
                s0 <- previous.best
                if (trace) 
                    cat("Starting from a previous solution\n")
            } else {
                init <- previous.best$points
                nc <- NCOL(init)
                if (nc > k)
                    init <- init[, 1:k, drop = FALSE]
                else   # nc < k
                    for (i in 1:(k-nc))
                        init <- cbind(init, runif(NROW(init), -0.1, 0.1))
                # evaluate isoMDS with stress
                s0 <- isoMDS(dist, init, k = k, maxit = 0)
                # zero 'tries': this was a new start
                if (inherits(previous.best, "metaMDS"))
                    previous.best$tries <- 0
                if (trace)
                    cat(gettextf("Starting from %d-dimensional solution\n", nc))
            }
        } else if (is.matrix(previous.best) || is.data.frame(previous.best)) {
            s0 <- isoMDS(dist, previous.best, k = k, maxit = 0)
            if (trace)
                cat("Starting from supplied configuration\n")
        } else { # an error!
            stop("'previous.best' of unknown kind")
        }
    } # No previous best:
    else s0 <- isoMDS(dist, k = k, trace = isotrace)
    if (trace) 
        cat("Run 0 stress", s0$stress, "\n")
    tries <- 0
    while(tries < trymax) {
        tries <- tries + 1
        stry <- isoMDS(dist, initMDS(dist, k = k), k = k, maxit = 200, 
                       tol = 1e-07, trace = isotrace)
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
    out <- list(points = s0$points, dims = k, stress = s0$stress, 
                data = attr(dist, "commname"),
                distance = attr(dist, "method"), converged = converged,
                tries = tries)
    out
}
