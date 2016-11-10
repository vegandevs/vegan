`intersetcor` <-
    function(object)
{
    if (!inherits(object, "cca"))
        stop("can be used only with objects inheriting from 'cca'")
    if (is.null(object$CCA))
        stop("can be used only with constrained ordination")
    wa <- object$CCA$wa
    X <- qr.X(object$CCA$QR)
    if (inherits(object, "rda"))
        cor(X, wa)
    else { # cca
        w <- object$rowsum
        cov <- crossprod(X, sqrt(w) * wa)
        sd <- outer(1/sqrt(colSums(X^2)), 1/sqrt(colSums(w * wa^2)))
        cov * sd
    }
}
