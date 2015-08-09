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
