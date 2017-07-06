## deprecated functions

## as.mlm: deprecated in 2.5-0

`as.mlm` <-
    function(x)
{
    .Deprecated("see ?hatvalues.cca for new alternatives")
    if (is.null(x$CCA))
        stop("'as.mlm' can be used only for constrained ordination")
    UseMethod("as.mlm")
}
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

`as.mlm.rda` <-
    function (x)
{
    X <- as.data.frame(qr.X(x$CCA$QR))
    WA <- x$CCA$wa
    lm(WA ~ . , data = X)
}

