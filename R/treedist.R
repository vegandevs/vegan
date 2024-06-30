`treedist` <-
    function(x, tree, relative = TRUE,  match.force = TRUE, ...)
{
    ## we cannot reconstruct tree with reversals from cophenetic
    tree <- as.hclust(tree)
    if (any(diff(tree$height) < -sqrt(.Machine$double.eps)))
        stop("tree with reversals cannot be handled")
    x <- as.matrix(x)
    n <- nrow(x)
    ABJ <- matrix(0, n , n)
    dmat <- as.matrix(cophenetic(tree))
    ## match names
    if (ncol(x) != ncol(dmat) || match.force) {
        if(!match.force)
            warning("dimensions do not match between 'x' and 'tree': matching by names")
        nm <- colnames(x)
        dmat <- dmat[nm, nm]
    }
    for(j in 1:n) {
        for (k in j:n) {
            jk <- x[j,] > 0 | x[k,] > 0
            if (sum(jk) > 1)
                ABJ[k, j] <- treeheight(update(tree, d = as.dist(dmat[jk, jk])))
        }
    }
    A <- diag(ABJ)
    AB <- as.dist(outer(A, A, "+"))
    ABJ <- as.dist(ABJ)
    out <- (2 * ABJ - AB)
    if (relative)
        out <- out/ABJ
    out[ABJ==0] <- 0
    attr(out, "method") <- if (relative) "treedist" else "raw treeedist"
    attr(out, "call") <- match.call()
    attr(out, "Labels") <- row.names(x)
    ## if (relative) theoretical maximum is 2, but that is only
    ## achieved when two zero-height trees (only one species) are
    ## combined into above zero-height tree (two species), and
    ## therefore we set here NA (but this can be reconsidered).
    attr(out, "maxdist") <- NA
    out
}
