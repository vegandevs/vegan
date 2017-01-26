`as.mlm.cca` <-
    function (x)
{
    w <- x$rowsum
    WA <- x$CCA$wa
    WA <- sweep(WA, 1, sqrt(w), "*")
    X <- qr.X(x$CCA$QR)
    lm(wa ~ ., data = as.data.frame(X))
}

