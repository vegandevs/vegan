`rda.default` <-
    function (X, Y = NULL, Z = NULL, scale = FALSE, ...)
{
    ## Protect against grave misuse: some people have used
    ## dissimilarities instead of data
    if (inherits(X, "dist") || NCOL(X) == NROW(X) &&
        isTRUE(all.equal(X, t(X))))
        stop("function cannot be used with (dis)similarities")
    X <- as.matrix(X)
    if (!is.null(Y)) {
        if (is.data.frame(Y) || is.factor(Y))
            Y <- model.matrix(~ ., as.data.frame(Y))[,-1,drop=FALSE]
        Y <- as.matrix(Y)
    }
    if (!is.null(Z)) {
        if (is.data.frame(Z) || is.factor(Z))
            Z <- model.matrix(~ ., as.data.frame(Z))[,-1,drop=FALSE]
        Z <- as.matrix(Z)
    }

    sol <- ordConstrained(X, Y, Z, arg = scale, method = "rda")

    call <- match.call()
    call[[1]] <- as.name("rda")
    sol$call <- call
    inertia <- if (scale) "correlations" else "variance"
    sol <- c(sol,
             list("inertia" = inertia))
    ## package klaR also has rda(): add a warning text that will be
    ## printed if vegan::rda object is displayed with klaR:::print.rda
    sol$regularization <- "this is a vegan::rda result object"
    class(sol) <- c("rda", "cca")
    sol
}
