`coef.cca` <-
    function (object, ...)
{
    if(is.null(object$CCA))
        stop("unconstrained models do not have coefficients")
    Q <- object$CCA$QR
    u <- object$CCA$u
    u <- sqrt(object$rowsum) * u
    qr.coef(Q, u)
}

