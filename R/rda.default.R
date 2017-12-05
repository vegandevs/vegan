`rda.default` <-
    function (X, Y = NULL, Z = NULL, scale = FALSE, ...)
{
    ## Protect against grave misuse: some people have used
    ## dissimilarities instead of data
    if (inherits(X, "dist") || NCOL(X) == NROW(X) &&
        isTRUE(all.equal(X, t(X))))
        stop("function cannot be used with (dis)similarities")
    X <- as.matrix(X)
    if (!is.null(Y))
        Y <- as.matrix(Y)
    if (!is.null(Z))
        Z <- as.matrix(Z)

    sol <- ordConstrained(X, Y, Z, arg = scale, method = "rda")

    call <- match.call()
    call[[1]] <- as.name("rda")
    sol$call <- call
    inertia <- if (scale) "correlations" else "variance"
    sol <- c(sol,
             list("inertia" = inertia))
    class(sol) <- c("rda", "cca")
    sol
}
