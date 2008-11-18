`metaMDSdist` <-
    function (comm, distance = "bray", autotransform = TRUE, noshare = 0.1, 
              trace = 1, commname, zerodist = "fail", distfun = vegdist, 
              ...) 
{
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
    if (any(dis <= 0)) {
        if (zerodist == "fail") 
            stop("Zero dissimilarities are not allowed")
        else if (zerodist == "add") {
            zero <- min(dis[dis > 0])/2
            dis[dis <= 0] <- zero
            warning("Zero dissimilarities changed into ", zero)
        }
    }
    maxdis <- max(dis)
    if (sum(tmp <- no.shared(comm))/length(dis) > noshare && noshare > 0) {
        if (trace) 
            cat("Using step-across dissimilarities:\n")
        dis <- stepacross(dis, trace = trace, ...)
    }
    if (length(unique(distconnected(tmp, trace = trace > 1))) > 1) 
        warning("Data are disconnected, results may be meaningless")
    attr(dis, "maxdis") <- maxdis
    attr(dis, "commname") <- commname
    attr(dis, "function") <- distname
    dis
}
