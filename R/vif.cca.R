`vif.cca` <-
    function(object) 
{
    Q <- object$CCA$QR
    out <- rep(NA, NCOL(Q$qr))
    names(out)[Q$pivot] <- colnames(Q$qr)
    rank <- Q$rank
    V <- chol2inv(Q$qr, size = rank)
    X <- qr.X(Q)[, Q$pivot[1:rank], drop=FALSE]
    Vi <- crossprod(X)
    v1 <- diag(V)
    v2 <- diag(Vi)
    out[Q$pivot[1:rank]] <- v1 * v2
    out
}

