`as.mlm.cca` <-
    function (x)
{
    w <- x$rowsum
    WA <- x$CCA$wa
    X <- qr.X(x$CCA$QR)
    ## shall use weighted regression: deweight X
    X <- (1/sqrt(w)) * X
    X <- as.data.frame(X)
    lm(WA ~ ., data = X, weights = w)
}

