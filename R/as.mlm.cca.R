`as.mlm.cca` <-
    function (x) 
{
    w <- x$rowsum
    wa <- x$CCA$wa
    wa <- sweep(wa, 1, sqrt(w), "*")
    X <- qr.X(x$CCA$QR)
    lm(wa ~ . - 1, data = as.data.frame(X))
}

