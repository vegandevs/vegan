"coef.cca" <-
function (object, ...) 
{
    if(is.null(object$CCA))
        stop("unconstrained models do not have coefficients")
    Q <- object$CCA$QR
    u <- object$CCA$u
    u <- sweep(u, 1, sqrt(object$rowsum), "*")
    qr.coef(Q, u)
}

