`betadisper` <-
    function(d, group, type = c("median","centroid"))
{
    ## inline function for spatial medians
    spatialMed <- function(vectors, group, pos) {
        axes <- seq_len(NCOL(vectors))
        spMedPos <- ordimedian(vectors, group, choices = axes[pos])
        spMedNeg <- ordimedian(vectors, group, choices = axes[!pos])
        return(cbind(spMedPos, spMedNeg))
    }
    ## Tolerance for zero Eigenvalues
    TOL <- 1e-7
    ## uses code from stats:::cmdscale by R Core Development Team
    if(!inherits(d, "dist"))
        stop("distances 'd' must be a 'dist' object")
    if(missing(type))
        type <- "median"
    type <- match.arg(type)
    ## checks for groups - need to be a factor for later
    if(!is.factor(group))
        group <- as.factor(group)
    n <- attr(d, "Size")
    x <- matrix(0, ncol = n, nrow = n)
    x[row(x) > col(x)] <- d^2
    ## site labels
    labs <- attr(d, "Labels")
    ## remove NAs in group
    if(any(gr.na <- is.na(group))) {
        group <- group[!gr.na]
        x <- x[!gr.na, !gr.na]
        ## update n otherwise C call crashes
        n <- n - sum(gr.na)
        ## update labels
        labs <- labs[!gr.na]
        warning("Missing observations due to 'group' removed.")
    }
    ## remove NA's in d
    if(any(x.na <- apply(x, 1, function(x) any(is.na(x))))) {
        x <- x[!x.na, !x.na]
        group <- group[!x.na]
        ## update n otherwise C call crashes
        n <- n - sum(x.na)
        ## update labels
        labs <- labs[!x.na]
        warning("Missing observations due to 'd' removed.")
    }
    x <- x + t(x)
    storage.mode(x) <- "double"
    .C("dblcen", x, as.integer(n), DUP = FALSE, PACKAGE="stats")
    e <- eigen(-x/2, symmetric = TRUE)
    vectors <- e$vectors
    eig <- e$values
    ## Remove zero eigenvalues
    eig <- eig[(want <- abs(eig/eig[1]) > TOL)]
    ## scale Eigenvectors
    vectors <- vectors[, want, drop = FALSE] %*% diag(sqrt(abs(eig)))
    ## store which are the positive eigenvalues
    pos <- eig > 0
    ## group centroids in PCoA space
    centroids <-
        switch(type,
               centroid = apply(vectors, 2, function(x) tapply(x, group, mean)),
               median = spatialMed(vectors, group, pos)
               )
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
        dist.pos <- sweep(vectors[, pos, drop=FALSE], 2, centroids[pos])
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
    colnames(vectors) <- names(eig) <- paste("PCoA", seq_along(eig), sep = "")
    if(is.matrix(centroids))
        colnames(centroids) <- names(eig)
    else
        names(centroids) <- names(eig)
    rownames(vectors) <- names(zij) <- labs
    retval <- list(eig = eig, vectors = vectors, distances = zij,
                   group = group, centroids = centroids, call = match.call())
    class(retval) <- "betadisper"
    attr(retval, "method") <- attr(d, "method")
    attr(retval, "type") <- type
    retval
}
