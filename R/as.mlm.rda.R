`as.mlm.rda` <-
    function (x) 
{
    X <- qr.X(x$CCA$QR)
    colnames(X) <- colnames(X)[x$CCA$QR$pivot]
    if (!is.null(x$na.action) && inherits(x$na.action, "exclude"))
        x$CCA$wa <- x$CCA$wa[-x$na.action,]
    lm(x$CCA$wa ~ . - 1, data = as.data.frame(X))
}

