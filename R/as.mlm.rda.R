`as.mlm.rda` <-
    function (x) 
{
    if (inherits(x, "dbrda"))
        stop("cannot be used with 'dbrda'")
    X <- qr.X(x$CCA$QR)
    lm(x$CCA$wa ~ . - 1, data = as.data.frame(X))
}

