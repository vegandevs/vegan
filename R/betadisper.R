`betadisper` <-
    function(d, group, type = c("centroid","median"))
{
    ## uses code from stats:::cmdscale by R Core Development Team
    if(!inherits(d, "dist"))
        stop("distances 'd' must be a 'dist' object")
    ## checks for groups - need to be a factor for later
    if(!is.factor(group))
        group <- as.factor(group)
    n <- attr(d, "Size")
    x <- matrix(0, ncol = n, nrow = n)
    x[row(x) > col(x)] <- d^2
    x <- x + t(x)
    .C(stats:::R_dblcen, as.double(x), as.integer(n), DUP = FALSE)
    e <- eigen(-x/2, symmetric = TRUE)
    vectors <- e$vectors
    eig <- e$values
    ## store which are the positive eigenvalues
    pos <- eig > 0
    ## scale Eigenvectors
    vectors <- vectors %*% diag(sqrt(abs(eig)))
    ## group centroids in PCoA space
    centroids <- apply(vectors, 2, function(x) tapply(x, group, mean))
    ## for each of the groups, calculate distance to centroid for
    ## observation in the group
    dist.pos <- vectors[, pos, drop=FALSE] - centroids[group, pos, drop=FALSE]
    dist.pos <- rowSums(dist.pos^2)
    if (any(!pos)) {
        dist.neg <- vectors[, !pos, drop=FALSE] -
            centroids[group, !pos, drop=FALSE]
        dist.neg <- rowSums(dist.neg^2)
    } else {
        dist.neg <- 0
    }
    ## zij are the distances of each point to its group centroid
    zij <- sqrt(abs(dist.pos - dist.neg))
    ## add in correct labels
    colnames(vectors) <- colnames(centroids) <- names(eig) <-
        paste("PCoA", 1:n, sep = "")
    rownames(vectors) <- names(zij) <- attr(d, "Labels")
    retval <- list(eig = eig, vectors = vectors, distances = zij,
                   group = group, centroids = centroids, call = match.call())
    class(retval) <- "betadisper"
    attr(retval, "method") <- attr(d, "method")
    retval
}
