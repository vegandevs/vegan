`coef.cca` <-
    function (object, norm = FALSE, ...)
{
    if(is.null(object$CCA))
        stop("unconstrained models do not have coefficients")
    Q <- object$CCA$QR
    u <- object$CCA$u
    ## if rank==0, the next would fail, but this kluge gives
    ## consistent results with coef.rda and vegan 2.4
    if (ncol(u))
        u <- sqrt(object$rowsum) * u
    ## scores.cca uses na.predict and may add missing NA rows to u,
    ## but Q has no missing cases
    if (nrow(Q$qr) < nrow(u) && inherits(object$na.action, "exclude"))
        u <- u[-object$na.action,, drop=FALSE]
    b <- qr.coef(Q, u)
    if (norm)
        b <- sqrt(colSums(qr.X(Q)^2)) * b
    b
}

