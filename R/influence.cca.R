### influence statistics for cca objects

## extract QR decomposition

`qr.cca` <-
    function(x, ...)
{
    if (is.null(r <- x$CCA$QR))
        stop("unconstrained model or rank zero: no 'qr' component")
    r
}

## residual degrees of freedom

`df.residual.cca` <-
    function(object, ...)
{
    nobs(object) - max(object$pCCA$rank, 0) - max(object$CCA$qrank, 0) - 1
}

## hat values need adjustment, because QR ignores Intercept

`hatvalues.cca` <-
    function(model, ...)
{
    hat <- rowSums(qr.Q(qr(model))^2) + weights(model)
    ## it can be be that hat > 1: the trick below is the same as used
    ## in stats::lm.influence()
    hat[hat > 1 - 10 * .Machine$double.eps] <- 1
    hat
}

`hatvalues.rda` <-
    function(model, ...)
{
    hat <- rowSums(qr.Q(qr(model))^2) + 1/nrow(qr(model)$qr)
    hat[hat > 1 - 10 * .Machine$double.eps] <- 1
    hat
}

## sigma gives the residual standard deviation. The only unambiguous
## sigma is the residual deviation for species, but for CANOCO like
## statistic this would be residual of WA/LC regression with little or
## no meaning.

`sigma.cca` <-
    function(object, type = c("response", "canoco"), ...)
{
    type <- match.arg(type)
    ## no response type in distance-based ordination
    if (inherits(object, c("dbrda", "capscale")) &&
        type == "response")
        stop("type = 'response' is not available in distance-based ordination")
    ## response: a vector of species (column) sigmata
    rdf <- df.residual(object)
    if (inherits(object, "rda"))
        adj <- nobs(object) - 1
    else
        adj <- 1
    if (type == "response") {
        sqrt(colSums(ordiYbar(object, "CA")^2 / rdf * adj))
    } else { # canoco has WA - LC regression
        sqrt((colSums(weights(object) * object$CCA$wa^2) - 1)/rdf)
    }
}

## rstandard and rstudent need sigma and have similar restrictions as
## sigma: it should be extractable and meaningful.

`rstandard.cca` <-
    function(model, type = c("response", "canoco"), ...)
{
    type <- match.arg(type)
    sd <- sigma(model, type = type)
    hat <- hatvalues(model)
    if (inherits(model, "rda"))
        adj <- sqrt(nobs(model) - 1)
    else
        adj <- 1
    w <- sqrt(weights(model))
    res <- switch(type,
                  "response" = ordiYbar(model, "CA") * adj,
                  "canoco" = w * (model$CCA$wa - model$CCA$u))
    res <- res / sqrt(1 - hat)
    res <- sweep(res, 2, sd, "/")
    res[is.infinite(res)] <- NaN
    attributes(res) <- list(dim = dim(res), dimnames = dimnames(res))
    res
}

## MASS book (4th ed), Chapter 6.3, p. 152 (e^star)

`rstudent.cca` <-
    function(model, type = c("response", "canoco"), ...)
{
    type <- match.arg(type)
    np <- df.residual(model)
    res <- rstandard(model, type = type)
    res <- res / sqrt(pmax.int(np-res^2, 0)/(np-1))
    res[is.infinite(res)] <- NaN
    attributes(res) <- list(dim = dim(res), dimnames = dimnames(res))
    res
}

## Cook's distance depends on meaningful sigma

`cooks.distance.cca` <-
    function(model, type = c("response", "canoco"), ...)
{
    hat <- hatvalues(model)
    p <- model$CCA$qrank
    d <- rstandard(model, type = type)^2 * hat / (1 - hat) / (p + 1)
    d[is.infinite(d)] <- NaN
    attributes(d) <- list(dim = dim(d), dimnames = dimnames(d))
    d
}

## residual sums of squares and products

`SSD.cca` <-
    function(object, type = "canoco", ...)
{
    type <- match.arg(type)
    w <- sqrt(weights(object))
    SSD <- crossprod(w * (object$CCA$wa - object$CCA$u))
    structure(list(SSD = SSD, call = object$call, df = df.residual(object)),
              class = "SSD")
}

## variances and covariances of coefficients. The sqrt(diag()) will be
## standard errors of regression coefficients, and for constrained
## ordination model m, the t-values of regression coefficients will be
## coef(m)/sqrt(diag(vcov(m))). The kind of relevant coefficient will
## be determined by the output of SSD.

`vcov.cca` <-
    function(object, type = "canoco", ...)
{
    type <- match.arg(type)
    QR <- qr(object)
    p <- 1L:QR$rank
    ## we do not give the (Intercept): it is neither in coef()
    cov.unscaled <- chol2inv(QR$qr[p, p])
    dimnames(cov.unscaled) <- list(colnames(QR$qr)[p], colnames(QR$qr)[p])
    ssd <- estVar(SSD(object, type = type))
    kronecker(ssd, cov.unscaled, make.dimnames = TRUE)
}
