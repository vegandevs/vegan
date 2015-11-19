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
            stop("Zero dissimilarities are not allowed")
        else if (zerodist == "add") {
            zero <- min(dis[dis > 0], na.rm = TRUE)/2
            dis[dis <= 0] <- zero
            if (trace)
                cat("Zero dissimilarities changed into ", zero,"\n")
        }
    }
    ## We actually used maxdis to decide whether index has a closed
    ## upper limit, but simple maximum does not give that info.
    ## Therefore we see if an arbitrary matrix with no shared species
    ## has distance = 1.
    maxdis <- abs(distfun(matrix(c(7,0,0,3), 2, 2),
                      method = distance, ...) - 1) < 1e-4
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
            warning("Data are disconnected, results may be meaningless")
    }
    attr(dis, "maxdis") <- maxdis
    attr(dis, "commname") <- commname
    attr(dis, "comm") <- comm
    attr(dis, "function") <- distname
    dis
}
