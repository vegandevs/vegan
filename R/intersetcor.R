`intersetcor` <-
    function(object)
{
    if (!inherits(object, "cca"))
        stop("can be used only with objects inheriting from 'cca'")
    if (is.null(object$CCA))
        stop("can be used only with constrained ordination")
    wa <- object$CCA$wa
    X <- qr.X(object$CCA$QR)
    ## remove conditions (partial terms)
    if (!is.null(object$pCCA)) {
        X <- X[, -seq_along(object$pCCA$envcentre), drop = FALSE]
        X <- qr.resid(object$pCCA$QR, X)
    }
    if (inherits(object, "rda"))
        cor(X, wa)
    else { # cca: weighted analysis, terms already weighted-centred
        wa <- sqrt(object$rowsum) * wa
        cov <- crossprod(X, wa)
        isd <- outer(1/sqrt(colSums(X^2)), 1/sqrt(colSums(wa^2)))
        cov * isd
    }
}
