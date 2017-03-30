### influence statistics for cca objects

## extract QR decomposition

`qr.cca` <-
    function(x, ...)
{
    if (is.null(r <- x$CCA$QR))
        stop("unconstrained model or rank zero: no 'qr' component")
    r
}

## hat values need adjustment, because QR ignores Intercept

`hatvalues.cca` <-
    function(model, ...)
{
    rowSums(qr.Q(model$CCA$QR)^2) + 1/nrow(model$CCA$QR$qr)
}
