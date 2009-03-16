`oecosimu` <-
    function(comm, nestfun, method, nsimul=99,
             burnin=0, thin=1, statistic = "statistic",
             control, ...)
{
    nestfun <- match.fun(nestfun)
    method <- match.arg(method, c("r00", "r0", "r1", "r2", "c0",
                                  "swap", "tswap", "backtrack",
                                  "quasiswap", "permat"))   # "permat" method added
    ## eveluate according to method
    if (method == "permat") {
        quant <- TRUE
        if (missing(control))
            control <- permat.control()
        pfull <- control$ptype == "full"
#        if (control$method %in% c("swap", "tswap", "abuswap")) {
#            if (thin != control$thin)
#                warning("'thin' and 'control$thin' not equal")
#            if (burnin != control$burnin)
#                warning("'burnin' and 'control$burnin' not equal")
#        }
    } else quant <- FALSE

    ## conditional on quant value
    if (!quant)
        comm <- ifelse(comm > 0, 1, 0)

    ind <- nestfun(comm, ...)
    if (is.list(ind))
        indstat <- ind[[statistic]]
    else
        indstat <- ind
    n <- length(indstat)
    simind <- matrix(0, nrow=n, ncol=nsimul)

    ## quant = FALSE
    if (!quant) {
        if (method %in% c("swap", "tswap")){
            checkbrd <- 1
            if (method == "tswap") {
                checkbrd <- sum(designdist(comm, "(J-A)*(J-B)", "binary"))
                M <- ncol(comm)
                N <- nrow(comm)
                checkbrd <- M*(M-1)*N*(N-1)/4/checkbrd
                thin <- round(thin*checkbrd)
            }
            attr(simind, "thin") <- thin
            attr(simind, "burnin") <- burnin
            x <- comm
            if (burnin > 0)
                x <- commsimulator(x, method= method, thin = round(checkbrd) * burnin)
            for(i in 1:nsimul) {
                x <- commsimulator(x, method = method, thin = thin)
                tmp <- nestfun(x, ...)
                if (is.list(tmp))
                    simind[,i] <- tmp[[statistic]]
                else
                    simind[,i] <- tmp
            }
        }
        else {
            for (i in 1:nsimul) {
                x <- commsimulator(comm, method=method)
                tmp <- nestfun(x,...)
                if (is.list(tmp))
                    simind[,i] <- tmp[[statistic]]
                else
                    simind[,i] <- tmp
            }
        }
    ## this is new addition for quantitative null model simulations
    ## quant = TRUE
    } else {
        ## permatfull
        if (pfull) {
            for (i in 1:nsimul) {
                x <- permatfull(comm, fixedmar=control$fixedmar,
                    shuffle=control$shuffle,
                    strata=control$strata,
                    mtype=control$mtype, times=1)
                tmp <- nestfun(x$perm[[1]])
                if (is.list(tmp)) {
                    simind[, i] <- tmp[[statistic]]
                } else simind[, i] <- tmp
            }
            attr(simind, "thin") <- NULL
            attr(simind, "burnin") <- NULL
        ## permatswap
        } else {
            if (control$method %in% c("swap", "tswap", "abuswap") && burnin > 0) {
                m <- permatswap(comm, method=control$method,
                    fixedmar=control$fixedmar,
                    shuffle=control$shuffle,
                    strata=control$strata,
                    mtype=control$mtype, times=1, 
                    burnin=burnin, thin=0)$perm[[1]]
            }  else m <- comm
            for (i in 1:nsimul) {
                x <- permatswap(m, method=control$method,
                    fixedmar=control$fixedmar,
                    shuffle=control$shuffle,
                    strata=control$strata,
                    mtype=control$mtype, times=1,
                    burnin=0, thin=thin)
                tmp <- nestfun(x$perm[[1]], ...)
                if (is.list(tmp))
                    simind[, i] <- tmp[[statistic]]
                else simind[, i] <- tmp
            }
            attr(simind, "thin") <- thin
            attr(simind, "burnin") <- burnin
        }
    }
    ## end of addition

    z <- (indstat - rowMeans(simind))/apply(simind, 1, sd)
    p <- 2*pmin(rowSums(indstat > simind), rowSums(indstat < simind))
    p <- (p + 1)/(nsimul + 1)

    ## ADDITION: if z is NA then it is not correct to calculate p values
    ## try e.g. oecosimu(dune, sum, "permat")
    if (any(is.na(z)))
        p[is.na(z)] <- NA
    ## collapse method with control$method
    if (method == "permat")
        method <- paste("permat", control$method, sep=".")

    if (is.null(names(indstat)))
        names(indstat) <- statistic
    if (!is.list(ind))
        ind <- list(statistic = ind)
    ind$oecosimu <- list(z = z, pval = p, simulated=simind, method=method,
                         statistic = indstat)
    class(ind) <- c("oecosimu", class(ind))
    ind
}
