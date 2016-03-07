`as.mlm.rda` <-
    function (x) 
{
    X <- qr.X(x$CCA$QR)
    lm(x$CCA$wa ~ . - 1, data = as.data.frame(X))
}

