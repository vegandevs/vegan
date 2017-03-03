`coef.rda` <-
    function (object, norm = FALSE, ...)
{
    if(is.null(object$CCA))
        stop("unconstrained models do not have coefficients")
    Q <- object$CCA$QR
    b <- qr.coef(Q, object$CCA$u)
    if (norm)
        b <- sqrt(colSums(qr.X(object$CCA$QR)^2)) * b
    b
}

