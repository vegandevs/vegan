### tolerance method for decorana
###
### After rescaling, all values should be 1

`tolerance.decorana` <-
    function(x, data, choices = 1:4, which = c("sites", "species"),
             useN2 = TRUE, ...)
{
    if (missing(data))
        stop("original response data must be given")
    which <- match.arg(which)
    ## Native decorana scaling (sites are WA of species) does not
    ## allow useN2 with species (but this can be done after scaling of
    ## results, and therefore the code below is ready for this).
    if (useN2 && which == "species")
        warning("useN2 is not implemented for species")
    EPS <- sqrt(.Machine$double.eps)
    ## transform data like decorana did
    if (x$iweigh)
        data <- downweight(data, x$fraction)
    if (!is.null(x$before))
        stop("before/after not yet implemented")
    ## see if data are plausible given decorana solution
    if (nrow(data) != nrow(x$rproj))
        stop("'data' have wrong row dimension")
    if (ncol(data) != nrow(x$cproj))
        stop("'data' have wrong col dimension")
    ## check the first eigenvalue
    ev1 <- svd(initCA(data), nv=0, nu=0)$d[1]^2
    ev0 <- if (x$ira) x$evals[1] else x$evals.decorana[1]
    if (!isTRUE(all.equal(ev0, ev1, check.attributes=FALSE)))
        stop("'data' are not plausible given 'decorana' result")
    ## preliminaries over: start working
    res <- switch(which,
                  "sites" = x$rproj,
                  "species" = x$cproj)
    tot <- switch(which,
                  "sites" = rowSums(data),
                  "species" = colSums(data))
    ## N2 is constant for axes. Here and elsewhere we still handle
    ## species, since support can be added later.
    if (useN2 && which != "species") {
        y <- switch(which,
                    "sites" = data,
                    "species" = t(data))
        y <- (y / rowSums(y))^2
        N2 <- 1 / rowSums(y, na.rm = TRUE) # 1/H
        N2[abs(N2 - 1) < EPS] <- 1
        N2scaling <- sqrt(pmax(1 - 1/N2, 0))
    }
    ## go over axes
    for(i in choices) {
        X <- data * outer(x$rproj[,i], x$cproj[,i], "-")^2
        X[X < 0] <- 0
        if (which == "species")
            X <- t(X)
        res[,i] <- sqrt(rowSums(X)/tot)
        if (useN2 && which == "sites")
            res[,i] <- res[,i] / N2scaling
    }
    res <- res[,choices, drop=FALSE]
    res[!is.finite(res) | res < EPS] <- 0
    class(res) <- c("tolerance.decorana", "tolerance.cca", "tolerance", "matrix")
    attr(res, "which") <- which
    attr(res, "scaling") <- "decorana"
    attr(res, "N2") <- if (useN2 && which != "species") N2 else NA
    res
}
