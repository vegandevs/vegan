`as.mlm.rda` <-
    function (x) 
{
    X <- qr.X(x$CCA$QR)
    colnames(X) <- colnames(X)[x$CCA$QR$pivot]
    lm(x$CCA$wa ~ . - 1, data = as.data.frame(X))
}

