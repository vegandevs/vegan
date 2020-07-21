`pcnm` <- function(dis, threshold, w, dist.ret = FALSE) {
    ## square matrix to dist
    if ((is.matrix(dis) || is.data.frame(dis)) &&
        isSymmetric(unname(as.matrix(dis))))
        dis <- as.dist(dis)
    if (!inherits(dis, "dist"))
        stop("'dis' does not appear to be distances")
    EPS <- sqrt(.Machine$double.eps)
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
    if (NCOL(res$vectors))
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
