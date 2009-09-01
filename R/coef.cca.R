"coef.cca" <-
function (object, ...) 
{
    k <- !is.na(object$rowsum)
    Q <- object$CCA$QR
    u <- object$CCA$u[k,]
    u <- sweep(u, 1, sqrt(object$rowsum[k]), "*")
    qr.coef(Q, u)
}

