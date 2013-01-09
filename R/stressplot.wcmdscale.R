### stressplot() methods for eigenvector ordinations wcmdscale, rda,
### cca, capscale

`stressplot.wcmdscale` <-
    function(object, k = 2, pch,  p.col = "blue", l.col = "red", lwd = 2, ...)
{
    ## NB: Ignores weights
    
    ## Get the ordination distances in k dimensions
    if (k > NCOL(object$points))
        stop("'k' cannot exceed the number of real dimensions")
    odis <- dist(object$points[,1:k, drop = FALSE])
    ## Reconstitute the original observed distances
    dis <- dist(object$points)
    if (!is.null(object$negaxes))
        dis <- sqrt(dis^2 - dist(object$negaxes)^2)
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
    ## Not yet done for pRDA
    if (!is.null(object$pCCA))
        stop("not implemented yet for partial RDA")
    ## Normalized scores to reconstruct data
    u <- cbind(object$CCA$u, object$CA$u)
    ev <- c(object$CCA$eig, object$CA$eig)
    ## normalizing constant
    nr <- NROW(u)
    const <- sqrt(ev * (nr-1))
    ## Distances
    dis <- dist(u %*% diag(const))
    odis <- dist(u[,1:k, drop=FALSE] %*% diag(const[1:k], nrow = k))
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
    ## NB: Ignores row weights!
    
    ## Not yet done for pCCA
    if (!is.null(object$pCCA))
        stop("not implemented yet for partial CCA")
    ## Normalized scores to reconstruct data
    u <- cbind(object$CCA$u, object$CA$u)
    sev <- sqrt(c(object$CCA$eig, object$CA$eig))
    ## Distances
    dis <- dist(u %*% diag(sev))
    odis <- dist(u[,1:k, drop=FALSE] %*% diag(sev[1:k], nrow = k))
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
    ## Not yet done for pRDA
    if (!is.null(object$pCCA))
        stop("not implemented yet for partial analysis")
    ## Normalized scores to reconstruct data
    u <- cbind(object$CCA$u, object$CA$u)
    ev <- c(object$CCA$eig, object$CA$eig)
    ## normalizing constant
    const <- sqrt(ev)
    ## Distances
    dis <- dist(u %*% diag(const))
    if (!is.null(object$CA$imaginary.u.eig))
        dis <- sqrt(dis^2 - dist(object$CA$imaginary.u.eig)^2)
    if (!is.null(object$ac))
        dis <- dis - object$ac
    odis <- dist(u[,1:k, drop=FALSE] %*% diag(const[1:k], nrow = k))
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
