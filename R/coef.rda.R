"coef.rda" <-
function (object, ...) 
{
    if(is.null(object$CCA))
        stop("unconstrained models do not have coefficients")
    Q <- object$CCA$QR
    qr.coef(Q, object$CCA$u)
}

