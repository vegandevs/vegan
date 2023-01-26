`metaMDSdist` <-
    function (comm, distance = "bray", autotransform = TRUE,
              noshare = TRUE, trace = 1, commname,
              zerodist = "ignore", distfun = vegdist, ...)
{
    ## metaMDSdist should get a raw data matrix, but if it gets a
    ## 'dist' object return that unchanged and quit silently.
    if (inherits(comm, "dist")  ||
        ((is.matrix(comm) || is.data.frame(comm)) &&
             isSymmetric(unname(as.matrix(comm)))))
        return(comm)
    distname <- deparse(substitute(distfun))
    distfun <- match.fun(distfun)
    zerodist <- match.arg(zerodist, c("fail", "add", "ignore"))
    formals(distfun) <- c(formals(distfun), alist(... = ))
    formals(stepacross) <- c(formals(stepacross), alist(... = ))
    if (missing(commname))
        commname <- deparse(substitute(comm))
    xam <- max(comm)
    if (autotransform && xam > 50) {
        comm <- sqrt(comm)
        commname <- paste("sqrt(", commname, ")", sep = "")
        if (trace)
            cat("Square root transformation\n")
    }
    if (autotransform && xam > 9) {
        comm <- wisconsin(comm)
        commname <- paste("wisconsin(", commname, ")", sep = "")
        if (trace)
            cat("Wisconsin double standardization\n")
    }
    dis <- distfun(comm, method = distance, ...)
    call <- attr(dis, "call")
    call[[1]] <- as.name(distname)
    attr(dis, "call") <- call
    if (zerodist != "ignore" && any(dis <= 0, na.rm = TRUE)) {
        if (zerodist == "fail")
            stop("zero-value dissimilarities are not allowed")
        else if (zerodist == "add") {
            zero <- min(dis[dis > 0], na.rm = TRUE)/2
            dis[dis <= 0] <- zero
            if (trace)
                cat("Zero dissimilarities changed into ", zero,"\n")
        }
    }
    ## We actually used maxdiss to decide whether index has a closed
    ## upper limit, but data maximum does not give that
    ## info. vegan::vegdist returns the info as an attribute of dis,
    ## but for other distance functions we guess constant maximum with
    ## arbitrary data matrix. This test is known to fail in some
    ## cases, but better so than assuming wrong maxdist: for instance,
    ## stats::dist(x, method="canberra") has maxdist ncol(x) --
    ## vegan::vegdist(x, "canberra") has maxdist 1, and we voluntarily
    ## fail here with stats::dist.
    if (is.null(attr(dis, "maxdist"))) {
        mat <- matrix(c(1,0,0, 1,0,0, 0,7,0, 0,3,0, 0,0,0.2,0,0,10),
                      nrow=3)
        dmat <- distfun(mat, method = distance, ...)
        if (sd(dmat) < sqrt(.Machine$double.eps) &&
            max(dis) - max(dmat) < sqrt(.Machine$double.eps))
        {
            maxdis <- max(dmat)
            attr(dis, "maxdist") <- maxdis
            message("assuming that theoretical maximum distance is ", maxdis)
        } else {
            attr(dis, "maxdist") <- NA
        }
    }
    ## sanity check of dissimilarities: either similarities or failed
    ## logic above
    maxdis <- attr(dis, "maxdist")
    if (!is.null(maxdis) && is.numeric(maxdis)) {
        if (max(dis) > maxdis + sqrt(.Machine$double.eps)) {
            warning("some dissimilarities exceed expected maximum ", maxdis)
            attr(dis, "maxdist") <- NA
        }
        if(maxdis < sqrt(.Machine$double.eps))
            warning("perhaps you have similarities instead of dissimilarities?")
    }
    if ((isTRUE(noshare) && any(tmp <- no.shared(comm))) ||
        (!is.logical(noshare) && noshare >= 0 &&
         sum(tmp <- no.shared(comm))/length(dis) > noshare)) {
        if (trace)
            cat("Using step-across dissimilarities:\n")
        rn <- range(dis[tmp], na.rm = TRUE)
        if (rn[2]/rn[1] > 1.01)
            warning("non-constant distances between points with nothing shared\n",
                    "  stepacross may be meaningless: consider argument 'noshare=0'")
        is.na(dis) <- tmp
        dis <- stepacross(dis, trace = trace, toolong=0, ...)
        if (length(unique(distconnected(tmp, trace = trace))) > 1)
            warning("data are disconnected, results may be meaningless")
    }
    attr(dis, "commname") <- commname
    attr(dis, "comm") <- comm
    attr(dis, "function") <- distname
    dis
}
