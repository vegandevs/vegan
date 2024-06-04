`isomapdist` <-
    function(dist, epsilon, k, path="shortest", fragmentedOK = FALSE, ...)
{
    EPS <- 1e-5
    op <- options(warn = 2)
    on.exit(options(op))
    if (!inherits(dist, "dist"))
        dist <- as.dist(dist)
    options(op)
    method <- attr(dist, "method")
    if (missing(epsilon) && missing(k))
        stop("either epsilon or k must be given")
    if (!missing(epsilon) && !missing(k))
        message("both epsilon and k given, using epsilon")
    if (!missing(epsilon))
        dist[dist >= epsilon-EPS] <- NA
    else {
        dist <- as.matrix(dist)
        diag(dist) <- NA
        is.na(dist) <- apply(dist, 2, function(x)
                             x > x[order(x, na.last=TRUE)[k]])
        dist <- pmax(as.dist(dist), as.dist(t(dist)), na.rm = TRUE)
    }
    fragm <- distconnected(dist, toolong=0,  trace=FALSE)
    take <- NULL
    if (length(unique(fragm)) > 1) {
        if (fragmentedOK) {
            warning("data are fragmented: taking the largest fragment")
            take <- fragm == as.numeric(names(which.max(table(fragm))))
            dist <- as.dist(as.matrix(dist)[take,take])
        } else {
            stop("data are fragmented")
        }
    }
    net <- which(!is.na(dist))
    attr(dist, "method") <- method
    dist <- stepacross(dist, path = path, toolong = 0, trace = FALSE)
    if (any(is.na(dist))) {
        ## never get here: we selected non-fragmented subset
        stop("dissimilarities contain NA and cannot be analysed: file a bug report")
    }
    if (missing(epsilon)) {
        attr(dist, "criterion") <-"k"
        attr(dist, "critval") <- k
    }
    else {
        attr(dist, "criterion") <- "epsilon"
        attr(dist, "critval") <- epsilon
    }
    attr(dist, "method") <- paste(attr(dist, "method"), "isomap")
    attr(dist, "maxdist") <- NA
    attr(dist, "net") <- net
    attr(dist, "take") <- take
    attr(dist, "call") <- match.call()
    dist
}
