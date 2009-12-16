`oecosimu` <-
    function(comm, nestfun, method, nsimul=99,
             burnin=0, thin=1, statistic = "statistic",
             alternative = c("two.sided", "less", "greater"),
             ...)
{
    alternative <- match.arg(alternative)
    nestfun <- match.fun(nestfun)
    if (!is.function(method)) {
        method <- match.arg(method, c("r00", "r0", "r1", "r2", "c0",
                                  "swap", "tswap", "backtrack", "quasiswap",
                                  "r2dtable"))
        if (method == "r2dtable") {
            nr <- rowSums(comm)
            nc <- colSums(comm)
            permfun <- function(z) r2dtable(1, nr, nc)[[1]]
        }
    } else {
        permfun <- match.fun(method)
        method <- "custom"
    }
    quant <- method %in% c("r2dtable", "custom")

    ind <- nestfun(comm, ...)
    if (is.list(ind))
        indstat <- ind[[statistic]]
    else
        indstat <- ind
    n <- length(indstat)
    simind <- matrix(0, nrow=n, ncol=nsimul)

    ## permutation for binary data
    if (!quant) {
        comm <- ifelse(comm > 0, 1, 0)
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
    ## permutation for count data
    } else {
        if (!all(dim(comm) == dim(permfun(comm))))
            stop("permutation function is not compatible with community matrix")
        ## sequential algorithms
        if (burnin > 0 || thin > 1) {
            if (burnin > 0) {
                m <- permfun(comm, burnin=burnin, thin=1)
            }  else m <- comm
            for (i in 1:nsimul) {
                tmp <- nestfun(permfun(m, burnin=0, thin=thin), ...)
                if (is.list(tmp))
                    simind[, i] <- tmp[[statistic]]
                else simind[, i] <- tmp
            }
            attr(simind, "thin") <- thin
            attr(simind, "burnin") <- burnin
        ## not sequential algorithms
        } else {
            for (i in 1:nsimul) {
                tmp <- nestfun(permfun(comm), ...)
                if (is.list(tmp)) {
                    simind[, i] <- tmp[[statistic]]
                } else simind[, i] <- tmp
            }
            attr(simind, "thin") <- NULL
            attr(simind, "burnin") <- NULL
        }
    }
    ## end of addition
    sd <- apply(simind, 1, sd)
    z <- (indstat - rowMeans(simind))/sd
    if (any(sd < sqrt(.Machine$double.eps)))
        z[sd < sqrt(.Machine$double.eps)] <- 0
    pless <- rowSums(indstat <= simind)
    pmore <- rowSums(indstat >= simind)
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
    if (method == "custom")
        attr(method, "permfun") <- permfun
    ind$oecosimu <- list(z = z, pval = p, simulated=simind, method=method,
                         statistic = indstat, alternative = alternative)
    class(ind) <- c("oecosimu", class(ind))
    ind
}

