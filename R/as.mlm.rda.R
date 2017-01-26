`as.mlm.rda` <-
    function (x)
{
    X <- qr.X(x$CCA$QR)
    WA <- x$CCA$wa
    lm(WA ~ . , data = as.data.frame(X))
}

