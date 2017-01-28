`intersetcor` <-
    function(object)
{
    if (!inherits(object, "cca"))
        stop("can be used only with objects inheriting from 'cca'")
    if (is.null(object$CCA))
        stop("can be used only with constrained ordination")
    lc <- object$CCA$u
    X <- qr.X(object$CCA$QR)
    if (inherits(object, "rda"))
        cor(X, lc)
    else { # cca: weighted analysis, terms already weighted-centred
        lc <- sqrt(object$rowsum) * lc
        cov <- crossprod(X, lc)
        isd <- outer(1/sqrt(colSums(X^2)), 1/sqrt(colSums(lc^2)))
        cov * isd
    }
}
