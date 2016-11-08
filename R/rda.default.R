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

    call <- match.call()
    call[[1]] <- as.name("rda")
    sol$call <- call
    inertia <- if (scale) "correlations" else "variance"
    sol <- c(sol,
             list(method = "rda", "inertia" = inertia))
    class(sol) <- c("rda", "cca")
    sol
}
