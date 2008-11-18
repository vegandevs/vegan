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
    if (!missing(previous.best) && !is.null(previous.best)) {
        s0 <- previous.best
        if (trace) 
            cat("Starting from a previous solution\n")
    }
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
    if (!missing(previous.best) && !is.null(previous.best$tries)) 
        tries <- tries + previous.best$tries
    out <- list(points = s0$points, dims = k, stress = s0$stress, 
                data = attr(dist, "commname"), distance = attr(dist, 
                                               "method"), converged = converged, tries = tries)
    out
}
