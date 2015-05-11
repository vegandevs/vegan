`pcnm` <- function(dis, threshold, w, dist.ret = FALSE) {
    if (!inherits(dis, "dist")) {
        dims <- dim(dis)
        if (length(unique(dims)) >1) {
            stop("'dis' does not appear to be a square distance matrix.")
        }
        dis <- as.dist(dis)
    }
    EPS <- sqrt(.Machine$double.eps)
    wa.old <- options(warn = -1)
    on.exit(options(wa.old))
    dis <- as.dist(dis)
    if (missing(threshold)) {
        threshold <- max(spantree(dis)$dist)
    }
    dis[dis > threshold] <- 4*threshold
    ## vegan:::wcmdscale is able to use weights which also means that
    ## 'k' need not be given, but all vecctors with >0 eigenvalues
    ## will be found
    mypcnm <- wcmdscale(dis, eig = TRUE, w=w)
    res <- list(vectors = mypcnm$points, values = mypcnm$eig,
                weights = mypcnm$weig)
    k <- ncol(mypcnm$points)
    res$vectors <- sweep(res$vectors, 2, sqrt(res$values[seq_len(k)]), "/")
    colnames(res$vectors) <- paste("PCNM", 1:k, sep="")
    res$threshold <- threshold
    if (dist.ret) {
        attr(dis, "method") <- paste(attr(dis, "method"), "pcnm")
        attr(dis, "threshold") <- threshold
        res$dist <- dis
    }
    class(res) <- "pcnm"
    res
}
