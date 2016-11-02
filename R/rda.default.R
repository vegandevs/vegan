`rda.default` <-
    function (X, Y = NULL, Z = NULL, scale = FALSE, ...)
{
    ## Protect against grave misuse: some people have used
    ## dissimilarities instead of data
    if (inherits(X, "dist") || NCOL(X) == NROW(X) &&
        isTRUE(all.equal(X, t(X))))
        stop("function cannot be used with (dis)similarities")
    X <- as.matrix(X)
    sol <- ordConstrained(X, Y, Z, scale = scale, method = "rda")
    ## back-scale Xbar/Fit and colsum to the scale of observations:
    ## this should be handled in support functions, but we test it
    ## here to see that the results are consistent with previus ones.
    if (!scale) {
        scl <- sqrt(nrow(X) - 1)
        if (!is.null(sol$pCCA))
            sol$pCCA$Fit <- scl * sol$pCCA$Fit
        if (!is.null(sol$CCA))
            sol$CCA$Xbar <- scl * sol$CCA$Xbar
        sol$CA$Xbar <- scl * sol$CA$Xbar
        sol$colsum <- scl * sol$colsum
    }
    call <- match.call()
    call[[1]] <- as.name("rda")
    sol$call <- call
    inertia <- if (scale) "correlations" else "variance"
    sol <- c(sol,
             list(method = "rda", "inertia" = inertia))
    class(sol) <- c("rda", "cca")
    sol
}
