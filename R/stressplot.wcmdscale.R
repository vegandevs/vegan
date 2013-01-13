### stressplot() methods for eigenvector ordinations wcmdscale, rda,
### cca, capscale

`stressplot.wcmdscale` <-
    function(object, k = 2, pch,  p.col = "blue", l.col = "red", lwd = 2, ...)
{
    ## Check that original distances can be reconstructed: this
    ## requires that all axes were calculated instead of 'k' first.
    hasdims <- NCOL(object$points)
    if (!is.null(object$negaxes))
        hasdims <- hasdims + NCOL(object$negaxes)
    if (hasdims < length(object$eig))
        stop("observed distances cannot be reconstructed: all axes were not calculated")
    ## Get the ordination distances in k dimensions
    if (k > NCOL(object$points))
        stop("'k' cannot exceed the number of real dimensions")
    w <- sqrt(object$weights)
    u <- diag(w) %*% object$points
    odis <- dist(u[,1:k, drop = FALSE])
    ## Reconstitute the original observed distances
    dis <- dist(u)
    if (!is.null(object$negaxes))
        dis <- sqrt(dis^2 - dist(diag(w) %*% object$negaxes)^2)
    ## additive constant is not implemented in wcmdscale (which
    ## returns 'ac = NA'), but the next statement would take care of
    ## that: we want to have the input distances as observed distances
    ## so that we need to subtract 'ac' here, although ordination
    ## distances 'odis' do not add up to 'dis' but to 'dis + ac'.
    if (!is.na(object$ac))
        dis <- dis - object$ac
    ##Plot
    if (missing(pch))
        if (length(dis) > 5000)
            pch <- "."
        else
            pch <- 1
    plot(dis, odis, pch = pch, col = p.col, xlab = "Observed Dissimilarity",
         ylab = "Ordination Distance", ...)
    abline(0, 1, col = l.col, lwd = lwd, ...)
    invisible(odis)
}

`stressplot.rda` <-
    function(object, k = 2, pch, p.col = "blue", l.col = "red", lwd = 2, ...)
{
    ## Normalized scores to reconstruct data
    u <- cbind(object$CCA$u, object$CA$u)
    ev <- c(object$CCA$eig, object$CA$eig)
    ## normalizing constant
    nr <- NROW(u)
    const <- sqrt(ev * (nr-1))
    u <- u %*% diag(const)
    ## Distances
    dis <- dist(cbind(u, object$pCCA$Fit))
    odis <- dist(cbind(u[,seq_len(k), drop=FALSE], object$pCCA$Fit))
    ## plot like above
        ## Plot
    if (missing(pch))
        if (length(dis) > 5000)
            pch <- "."
        else
            pch <- 1
    plot(dis, odis, pch = pch, col = p.col, xlab = "Observed Dissimilarity",
         ylab = "Ordination Distance", ...)
    abline(0, 1, col = l.col, lwd = lwd, ...)
    invisible(odis)
}

`stressplot.cca` <-
    function(object, k = 2, pch, p.col = "blue", l.col = "red", lwd = 2, ...)
{
    ## Normalized scores to reconstruct data
    u <- cbind(object$CCA$u, object$CA$u)
    sev <- sqrt(c(object$CCA$eig, object$CA$eig))
    w <- sqrt(object$rowsum)
    u <- diag(w) %*% u %*% diag(sev)
    ## Distances
    dis <- dist(cbind(u, object$pCCA$Fit))
    odis <- dist(cbind(u[,seq_len(k), drop = FALSE], object$pCCA$Fit))
    ##odis <- dist(sweep(Xbar, 2, sqrt(object$colsum), "*"))
    ## plot like above
        ## Plot
    if (missing(pch))
        if (length(dis) > 5000)
            pch <- "."
        else
            pch <- 1
    plot(dis, odis, pch = pch, col = p.col, xlab = "Observed Dissimilarity",
         ylab = "Ordination Distance", ...)
    abline(0, 1, col = l.col, lwd = lwd, ...)
    invisible(odis)
}

`stressplot.capscale` <-
    function(object, k = 2, pch, p.col = "blue", l.col = "red", lwd = 2, ...)
{
    ## Scores to reconstruct data
    u <- cbind(object$CCA$u, object$CA$u)
    ev <- c(object$CCA$eig, object$CA$eig)
    u <- u %*% diag(sqrt(ev))
    if (!is.null(object$pCCA))
        pFit <- object$pCCA$Fit/sqrt(nrow(object$pCCA$Fit) - 1)
    else
        pFit <- NULL
    ## Distances
    dis <- dist(cbind(u, pFit))
    if (!is.null(object$CA$imaginary.u.eig))
        dis <- sqrt(dis^2 - dist(object$CA$imaginary.u.eig)^2)
    if (!is.null(object$ac))
        dis <- dis - object$ac
    odis <- dist(cbind(u[,seq_len(k), drop=FALSE], pFit))
    ## plot like above
        ## Plot
    if (missing(pch))
        if (length(dis) > 5000)
            pch <- "."
        else
            pch <- 1
    plot(dis, odis, pch = pch, col = p.col, xlab = "Observed Dissimilarity",
         ylab = "Ordination Distance", ...)
    abline(0, 1, col = l.col, lwd = lwd, ...)
    invisible(odis)
}

## Standard R PCA functions

`stressplot.prcomp` <-
    function(object, k = 2, pch, p.col = "blue", l.col = "red", lwd = 2, ...)
{
    dis <- dist(object$x)
    odis <- dist(object$x[, 1:k, drop = FALSE])
    if (missing(pch))
        if (length(dis) > 5000)
            pch <- "."
        else
            pch <- 1
    plot(dis, odis, pch = pch, col = p.col, xlab = "Observed Dissimilarity",
         ylab = "Ordination Distance", ...)
    abline(0, 1, col = l.col, lwd = lwd, ...)
    invisible(odis)    
}

`stressplot.princomp` <-
    function(object, k = 2, pch, p.col = "blue", l.col = "red", lwd = 2, ...)
{
    dis <- dist(object$scores)
    odis <- dist(object$scores[, 1:k, drop = FALSE])
    if (missing(pch))
        if (length(dis) > 5000)
            pch <- "."
        else
            pch <- 1
    plot(dis, odis, pch = pch, col = p.col, xlab = "Observed Dissimilarity",
         ylab = "Ordination Distance", ...)
    abline(0, 1, col = l.col, lwd = lwd, ...)
    invisible(odis)    
}
