`as.mlm.cca` <-
    function (x)
{
    w <- x$rowsum
    WA <- x$CCA$wa
    WA <- sweep(WA, 1, sqrt(w), "*")
    X <- qr.X(x$CCA$QR)
    lm(WA ~ ., data = as.data.frame(X))
}

