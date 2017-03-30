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
    rowSums(qr.Q(qr(model))^2) + 1/nrow(qr(model)$qr)
}

## sigma gives the residual standard deviation. The only unambiguous
## sigma is the residual deviation for species, but for CANOCO like
## statistic this would be residual of WA/LC regression with little or
## no meaning.

`sigma.cca` <-
    function(object, ...)
{
    ## a vector of species (column) sigmata
    rdf <- nobs(object) - m$CCA$qrank - 1
    colSums(ordiYbar(object, "CA")^2/rdf)
}
