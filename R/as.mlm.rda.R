`as.mlm.rda` <-
    function (x) 
{
    X <- qr.X(x$CCA$QR)
    ## We hope that column names will be fixed in R 2.12.1 (but
    ## perhaps in vain).
    if (getRversion() <= "2.12.0")
        colnames(X)[x$CCA$QR$pivot] <- colnames(X)
    lm(x$CCA$wa ~ . - 1, data = as.data.frame(X))
}

