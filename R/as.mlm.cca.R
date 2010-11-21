`as.mlm.cca` <-
    function (x) 
{
    w <- x$rowsum
    wa <- x$CCA$wa
    wa <- sweep(wa, 1, sqrt(w), "*")
    X <- qr.X(x$CCA$QR)
    ## qr.X gives wrong column names now, and they are fixed here
    ## (hopefully fixed in 2.12.1, but that's only a wish).
    if (getRversion() <= "2.12.0")
        colnames(X)[x$CCA$QR$pivot] <- colnames(X)
    lm(wa ~ . - 1, data = as.data.frame(X))
}

