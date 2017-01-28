`as.mlm.rda` <-
    function (x)
{
    X <- as.data.frame(qr.X(x$CCA$QR))
    WA <- x$CCA$wa
    lm(WA ~ . , data = X)
}

