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
    v <- cbind(object$CCA$v, object$CA$v)
    ev <- c(object$CCA$eig, object$CA$eig)
    ## normalizing constant
    nr <- NROW(u)
    const <- sqrt(ev * (nr-1))
    u <- u %*% diag(const)
    ## Distances
    Xbar <- u %*% t(v)
    Xbark <- u[, seq_len(k), drop = FALSE] %*% t(v[, seq_len(k), drop = FALSE])
    if (!is.null(object$pCCA)) {
        Xbar <- Xbar + object$pCCA$Fit
        Xbark <- Xbark + object$pCCA$Fit
    }
    dis <- dist(Xbar)
    odis <- dist(Xbark)
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
    v <- cbind(object$CCA$v, object$CA$v)
    v <- diag(sqrt(object$colsum)) %*% v
    ## Distances
    Xbar <- u %*% t(v)
    Xbark <- u[,seq_len(k), drop = FALSE] %*% t(v[,seq_len(k), drop = FALSE])
    if (!is.null(object$pCCA)) {
        Xbar <- Xbar + object$pCCA$Fit
        Xbark <- Xbark + object$pCCA$Fit
    }
    dis <- dist(Xbar)
    odis <- dist(Xbark)
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
    if (object$adjust == 1)
        const <- sqrt(NROW(u) - 1)
    else
        const <- 1
    u <- u %*% diag(sqrt(ev) * const)
    ## Constrained ordination needs also scores 'v' to reconstruct
    ## 'data', but these are not returned by capscale() which replaces
    ## original 'v' with weighted sums of 'comm' data.
    if (!is.null(object$CCA)) 
        v <- svd(object$CCA$Xbar - object$CA$Xbar, nu = 0, nv = object$CCA$qrank)$v
    else
        v <- NULL
    if (!is.null(object$CA))
        v <- cbind(v, svd(object$CA$Xbar, nu = 0, nv = object$CA$rank)$v)
    ## Reconstruct Xbar and Xbark
    Xbar <- u %*% t(v)
    Xbark <- u[,seq_len(k), drop = FALSE] %*% t(v[,seq_len(k), drop = FALSE])
    if (!is.null(object$pCCA)) {
        pFit <- object$pCCA$Fit/object$adjust
        Xbar <- Xbar + pFit
        Xbark <- Xbark + pFit
    }
    ## Distances
    dis <- dist(Xbar)
    odis <- dist(Xbark)
    if (!is.null(object$CA$imaginary.u.eig))
        dis <- sqrt(dis^2 - dist(object$CA$imaginary.u.eig)^2)
    if (!is.null(object$ac))
        dis <- dis - object$ac
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
