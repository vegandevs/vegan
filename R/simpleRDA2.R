### An internal function used in varpart(): Returns only the raw
### Rsquare and the rank of constraints in RDA.
`simpleRDA2` <-
    function (Y, X, SS.Y, ...)
{
    Q <- qr(X, tol=1e-6)
    Yfit.X <- qr.fitted(Q, Y)
    SS <- sum(Yfit.X^2)
    if (missing(SS.Y)) SS.Y <- sum(Y^2)
    Rsquare <- SS/SS.Y
    list(Rsquare = Rsquare, m = Q$rank)
}

### Analogous function, but the input must be Gower double-centred
### dissimilarities 'G = -GowerDblcen(as.matrix(dist(Y)^2))/2'. The
### math is based on McArdle & Anderson, Ecology 82: 290-297 (2001).
`simpleDBRDA` <-
    function(G, X, SS.G, ...)
{
    Q <- qr(X, tol=1e-6)
    Yfit.X <- qr.fitted(Q, G)
    SS <- sum(diag(Yfit.X))
    if (missing(SS.G)) SS.G <- sum(diag(G))
    Rsquare <- SS/SS.G
    list(Rsquare = Rsquare, m = Q$rank)
}

### Redesign Legendre's dbRDA.D function using QR decomposition. d is
### a "dist" object and X is a model matrix

`dbRDA2` <-
    function(d, X)
{
    ## basic manipulation
    G <- -GowerDblcen(as.matrix(d^2))/2
    X <- scale(X)
    ## QR decomposition
    Q <- qr(X, tol = 1e-6)
    ## Collect goodness-of-fit statistics
    SS <- sum(diag(G))
    SSfit <- sum(diag(qr.fitted(Q, G)))
    ## H is the hat matrix, and constrained ordination of LC scores is
    ## the eigen solution of HGH
    H <- tcrossprod(qr.Q(Q))
    HGH <- H %*% G %*% H
    e <- eigen(HGH, symmetric = TRUE)
    ## output using similar names as dbRDA.D. pos will select only
    ## above-zero eigenvalues and corresponding eigenvectors
    pos <- e$values > sqrt(.Machine$double.eps)
    sol <- list(SS.total = SS, SS.fit = SSfit, RDA.values = e$values[pos],
                RDA.coord = e$vectors[,pos] %*% diag(sqrt(e$values[pos])))
    sol
}
