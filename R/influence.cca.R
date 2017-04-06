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
     rowSums(qr.Q(qr(model))^2) + weights(model)
 }

`hatvalues.rda` <-
    function(model, ...)
{
    rowSums(qr.Q(qr(model))^2) + 1/nrow(qr(model)$qr)
}

## sigma gives the residual standard deviation. The only unambiguous
## sigma is the residual deviation for species, but for CANOCO like
## statistic this would be residual of WA/LC regression with little or
## no meaning.

`sigma.cca` <-
    function(object, type = c("response", "canoco"), ...)
{
    type <- match.arg(type)
    ## response: a vector of species (column) sigmata
     if (type == "response") {
        sqrt(colSums(ordiYbar(object, "CA")^2))
    } else { # canoco has WA - LC regression
        sqrt(colSums(weights(object) * object$CCA$wa^2) - 1)
    }
}

`sigma.rda` <-
    function(object, type = c("response", "canoco"), ...)
{
    type <- match.arg(type)
    ## response: a vector of species (column) sigmata
    rdf <- nobs(object) - object$CCA$qrank - 1
    if (type == "response") {
        sqrt(colSums(ordiYbar(object, "CA")^2/rdf))
    } else { # canoco has WA - LC regression
        sqrt((colSums(object$CCA$wa^2) - 1)/rdf)
    }
}

## rstandard and rstudent need sigma and have similar restrictions as
## sigma: it should be extractable and meaningful.

`rstandard.rda` <-
    function(model, type = c("response", "canoco"), ...)
{
    type <- match.arg(type)
    sd <- sigma(model, type = type)
    hat <- hatvalues(model)
    res <- switch(type,
                  "response" = ordiYbar(model, "CA"),
                  "canoco" = model$CCA$wa - model$CCA$u)
    res <- res / sqrt(1 - hat)
    res <- sweep(res, 2, sd, "/")
    res
}

## MASS book (4th ed), Chapter 6.3, p. 152 (e^star)

`rstudent.rda` <-
    function(model, type = c("response", "canoco"), ...)
{
    type <- match.arg(type)
    np <- nobs(model) - model$CCA$qrank - 1 # -1: Intercept
    res <- rstandard(model, type = type)
    res / sqrt((np-res^2)/(np-1))
}

## Cook's distance depends on meaningful sigma

`cooks.distance.rda` <-
    function(model, type = c("response", "canoco"), ...)
 {
     hat <- hatvalues(model)
     p <- model$CCA$qrank
     rstandard(model, type = type)^2 * hat / (1 - hat) / (p + 1)
 }
