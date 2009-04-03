`betadisper` <-
    function(d, group, type = c("centroid","median"), tol = 1e-07)
{
    ## uses code from stats:::cmdscale by R Core Development Team
    if(!inherits(d, "dist"))
        stop("distances 'd' must be a 'dist' object")
    type <- match.arg(type)
    if(type == "median")
        .NotYetUsed('type = "median"')
    ## checks for groups - need to be a factor for later
    if(!is.factor(group))
        group <- as.factor(group)
    n <- attr(d, "Size")
    x <- matrix(0, ncol = n, nrow = n)
    x[row(x) > col(x)] <- d^2
    x <- x + t(x)
    storage.mode(x) <- "double"
    .C("dblcen", x, as.integer(n), DUP = FALSE, PACKAGE="stats")
    e <- eigen(-x/2, symmetric = TRUE)
    vectors <- e$vectors
    eig <- e$values
    ## check d is Euclidean
    w0 <- eig[n] / eig[1]
    if(w0 > -tol)
        r <- sum(eig > (eig[1] * tol))
    else
        r <- length(eig)
    ## truncate eig if d is Euclidean
    eig <- eig[(rs <- seq_len(r))]
    ## scale Eigenvectors
    vectors <- vectors[, rs, drop = FALSE] %*% diag(sqrt(abs(eig)))
    ## store which are the positive eigenvalues
    pos <- eig > 0
    ## group centroids in PCoA space
    centroids <- apply(vectors, 2, function(x) tapply(x, group, mean))
    ## for each of the groups, calculate distance to centroid for
    ## observation in the group
    if(is.matrix(centroids)) {
        dist.pos <- vectors[, pos, drop=FALSE] -
            centroids[group, pos, drop=FALSE]
        dist.pos <- rowSums(dist.pos^2)
        if (any(!pos)) {
            dist.neg <- vectors[, !pos, drop=FALSE] -
                centroids[group, !pos, drop=FALSE]
            dist.neg <- rowSums(dist.neg^2)
        } else {
            dist.neg <- 0
        }
    } else {
        dist.pos <- sweep(vectors[, pos, drop=FALSE], 2,
                          centroids[pos])
        dist.pos <- rowSums(dist.pos^2)
        if (any(!pos)) {
            dist.neg <- sweep(vectors[, !pos, drop=FALSE], 2,
                              centroids[!pos])
            dist.neg <- rowSums(dist.neg^2)
        } else {
            dist.neg <- 0
        }
    }
    ## zij are the distances of each point to its group centroid
    zij <- sqrt(abs(dist.pos - dist.neg))
    ## add in correct labels
    colnames(vectors) <- names(eig) <- paste("PCoA", rs, sep = "")
    if(is.matrix(centroids))
        colnames(centroids) <- names(eig)
    else
        names(centroids) <- names(eig)
    rownames(vectors) <- names(zij) <- attr(d, "Labels")
    retval <- list(eig = eig, vectors = vectors, distances = zij,
                   group = group, centroids = centroids, call = match.call())
    class(retval) <- "betadisper"
    attr(retval, "method") <- attr(d, "method")
    retval
}
