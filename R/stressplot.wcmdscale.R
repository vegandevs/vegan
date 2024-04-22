### stressplot() methods for eigenvector ordinations wcmdscale, rda,
### cca, capscale, dbrda

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
        warning(gettextf("max allowed rank is k = %d", NCOL(object$points)))
    k <- min(NCOL(object$points), k)
    w <- sqrt(object$weights)
    u <- diag(w) %*% object$points
    odis <- dist(u[,1:k, drop = FALSE])
    ## Reconstitute the original observed distances
    dis <- dist(u)
    if (!is.null(object$negaxes))
        dis <- sqrt(dis^2 - dist(diag(w) %*% object$negaxes)^2)
    ## Remove additive constant to get original dissimilarities
    if (!is.na(object$ac)) {
        if (object$add == "lingoes")
            dis <- sqrt(dis^2 - 2 * object$ac)
        else if (object$add == "cailliez")
            dis <- dis - object$ac
        else
            stop("unknown Euclidifying adjustment: no idea what to do")
    }
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
    ## check that k does not exceed rank
    if (k > length(ev)) {
        warning(gettextf("max allowed rank is k = %d", length(ev)))
        k <- min(k, length(ev))
    }
    ## normalizing constant
    nr <- NROW(u)
    const <- sqrt(ev * (nr-1))
    u <- u %*% diag(const, length(const))
    ## Distances
    Xbar <- u %*% t(v)
    Xbark <- u[, seq_len(k), drop = FALSE] %*% t(v[, seq_len(k), drop = FALSE])
    if (!is.null(object$pCCA)) {
        pFit <- ordiYbar(object, "pCCA") * sqrt(nr-1)
        Xbar <- Xbar + pFit
        Xbark <- Xbark + pFit
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
    u <- diag(w) %*% u %*% diag(sev, length(sev))
    v <- cbind(object$CCA$v, object$CA$v)
    v <- diag(sqrt(object$colsum)) %*% v
    ## check that k <= rank
    if (k > length(sev)) {
        warning(gettextf("max allowed rank is k = %d", length(sev)))
        k <- min(k, length(sev))
    }
    ## Distances
    Xbar <- u %*% t(v)
    Xbark <- u[,seq_len(k), drop = FALSE] %*% t(v[,seq_len(k), drop = FALSE])
    if (!is.null(object$pCCA)) {
        pFit <- ordiYbar(object, "pCCA")
        Xbar <- Xbar + pFit
        Xbark <- Xbark + pFit
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
    ## check rank
    if (k > NCOL(u))
        warning(gettextf("max allowed rank is k = %d", ncol(u)))
    k <- min(k, ncol(u))
    ev <- c(object$CCA$eig, object$CA$eig)
    u <- u %*% diag(sqrt(ev) * object$adjust, length(ev))
    ## Constrained ordination needs also scores 'v' to reconstruct
    ## 'data', but these are not returned by capscale() which replaces
    ## original 'v' with weighted sums of 'comm' data.
    if (!is.null(object$CCA))
        v <- svd(ordiYbar(object, "CCA"), nu = 0, nv = object$CCA$qrank)$v
    else
        v <- NULL
    if (!is.null(object$CA))
        v <- cbind(v, svd(ordiYbar(object, "CA"), nu = 0,
                          nv = object$CA$rank)$v)
    ## Reconstruct Xbar and Xbark
    Xbar <- u %*% t(v)
    Xbark <- u[,seq_len(k), drop = FALSE] %*% t(v[,seq_len(k), drop = FALSE])
    if (!is.null(object$pCCA)) {
        pFit <- ordiYbar(object, "pCCA")
        Xbar <- Xbar + pFit
        Xbark <- Xbark + pFit
    }
    ## Distances
    dis <- dist(Xbar)
    odis <- dist(Xbark)
    if (!is.null(object$CA$imaginary.u.eig)) {
        dis <- dis^2 - dist(object$CA$imaginary.u.eig)^2
        if (all(dis > -sqrt(.Machine$double.eps)))
            dis <- sqrt(pmax.int(dis, 0))
        else # neg dis will be NaN with a warning
            dis <- sqrt(dis)
    }
    ## Remove additive constant to get original dissimilarities
    if (!is.null(object$ac)) {
        if (object$add == "lingoes")
            dis <- sqrt(dis^2 - 2 * object$ac)
        else if (object$add == "cailliez")
            dis <- dis - object$ac
        else
            stop("unknown Euclidifying adjustment: no idea what to do")
    }
    ## undo internal sqrt.dist
    if (object$sqrt.dist)
        dis <- dis^2
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

### dbrda() returns only row scores 'u' (LC scores for constraints,
### site scores for unconstrained part), and these can be used to
### reconstitute dissimilarities only in unconstrained ordination or
### for constrained component. stressplot across component would need
### reconstruction of data, but dbrda components are not additive, or
### ordiYbar(x, "initial") is *not* ordiYbar(x, "pCCA") + ordiYbar(x,
### "CCA") + ordiYbar(x, "CA").

`stressplot.dbrda` <-
    function(object, k = 2, pch, p.col = "blue", l.col = "red", lwd = 2, ...)
{
    if (!is.null(object$pCCA))
        stop("partial models cannot be analysed")
    ## Reconstructed zero distances can be tiny (negative) non-zero
    ## values, and we zap them to zero
    ZAP <- sqrt(.Machine$double.eps)
    ## Reconstruct original distances from Gower 'G'
    dis <- ordiYbar(object, "initial")
    dia <- diag(dis)
    dis <- -2 * dis + outer(dia, dia, "+")
    dis[abs(dis) < ZAP] <- 0
    dis <- sqrt(as.dist(dis)) * object$adjust
    ## Remove additive constant to get original dissimilarities
    if (!is.null(object$ac)) {
        if (object$add == "lingoes")
            dis <- sqrt(dis^2 - 2 * object$ac)
        else if (object$add == "cailliez")
            dis <- dis - object$ac
        else
            stop("unknown Euclidifying adjustment: no idea what to do")
    }
    ## undo internal sqrt.dist
    if (object$sqrt.dist)
        dis <- dis^2
    ## Approximate dissimilarities from real components with positive
    ## eigenvalues.
    if (!is.null(object$CCA)) {
        U <- object$CCA$u
        kmax <- object$CCA$poseig
        eig <- object$CCA$eig[seq_len(kmax)]
    } else {
        U <- object$CA$u
        kmax <- object$CA$poseig
        eig <- object$CA$eig[seq_len(kmax)]
    }
    U <- U %*% diag(sqrt(eig), nrow = kmax)
    if (k > kmax) {
        warning(gettextf("max allowed rank is k = %d", kmax))
        k <- kmax
    }
    if (k > 0) {
        odis <- dist(U[, seq_len(k), drop = FALSE]) * object$adjust
    } else {
        odis <- dist(matrix(0, nrow(U)))
    }
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
