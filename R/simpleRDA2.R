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
    R2adj <- RsquareAdj(Rsquare, nrow(Y), Q$rank)
    list(Rsquare = Rsquare, RsquareAdj = R2adj, m = Q$rank)
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
    R2adj <- RsquareAdj(Rsquare, nrow(G), Q$rank)
    list(Rsquare = Rsquare, RsquareAdj = R2adj, m = Q$rank)
}

### Analogous function for CCA. We initialize data with weighted
### double standaradization, and centre constraints X by row
### weights. The approximation of weighted R-square is found via
### permutations in permat (which must be given).

`simpleCCA` <-
    function(Y, X, SS.Y, permat, ...)
{
    Y <- initCA(Y)
    if(missing(SS.Y)) SS.Y <- sum(Y^2)
    w <- attr(Y, "RW")
    X <- .Call(do_wcentre, X, w, PACKAGE = "vegan")
    Q <- qr(X, tol=1e-6)
    Yfit.X <- qr.fitted(Q, Y)
    SS <- sum(Yfit.X^2)
    Rsquare <- SS/SS.Y
    ## permutation to estimate adjusted R2
    meanSS <- mean(sapply(seq_len(nrow(permat)),
                          function(i) sum(qr.fitted(Q, Y[permat[i,],])^2)))
    R2adj <- 1 - ((1 - Rsquare) / (1 - meanSS/SS.Y))
    list(Rsquare = Rsquare, RsquareAdj = R2adj, m = Q$rank)
}
