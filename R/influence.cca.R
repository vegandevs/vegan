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

`sigma.rda` <-
    function(object, ...)
{
    ## a vector of species (column) sigmata
    rdf <- nobs(object) - object$CCA$qrank - 1
    sqrt(colSums(ordiYbar(object, "CA")^2/rdf))
}

## rstandard and rstudent need sigma and have similar restrictions as
## sigma: it should be extractable and meaningful.

`rstandard.rda` <-
    function(model, ...)
{
    sd <- sigma(model)
    hat <- hatvalues(model)
    ## implement for working residuals: hardly interesting
    res <- ordiYbar(model, "CA")
    res <- res / sqrt(1 - hat)
    res <- sweep(res, 2, sd, "/")
    res
}

## MASS book e^star

`rstudent.rda` <-
    function(model, ...)
{
    np <- nobs(model) - model$CCA$qrank - 1 # -1: Intercept
    res <- rstandard(model)
    res / sqrt((np-res^2)/(np-1))
}

## Cook's distance depends on meaningful sigma

`cooks.distance.rda` <-
    function(model, ...)
 {
     hat <- hatvalues(model)
     p <- model$CCA$qrank
     res <- residuals(model, type="working") * sqrt(nobs(model)-p-1)
     ((res/(sigma(model) * (1-hat)))^2 * hat) / p
 }
