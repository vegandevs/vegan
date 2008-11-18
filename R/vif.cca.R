"vif.cca" <-
function (object) 
{
    Q <- object$CCA$QR
    V <- chol2inv(Q$qr)
    X <- qr.X(Q)
    Vi <- crossprod(X)
    nam <- colnames(Vi)
    v1 <- diag(V)
    v2 <- diag(Vi)
    structure(v1 * v2, names = nam)
}

