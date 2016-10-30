`cca.default` <-
    function (X, Y = NULL, Z = NULL, ...)
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
    if (any(rowSums(X) <= 0))
        stop("All row sums must be >0 in the community data matrix")
    if (any(tmp <- colSums(X) <= 0)) {
        exclude.spec <- seq(along=tmp)[tmp]
        names(exclude.spec) <- colnames(X)[tmp]
        class(exclude.spec) <- "exclude"
        X <- X[, !tmp, drop = FALSE]
    }
    sol <- ordConstrained(X, Y, Z, method = "cca")
    if (exists("exclude.spec")) {
        attr(sol$CA$v, "na.action") <- exclude.spec
    }
    call <- match.call()
    call[[1]] <- as.name("cca")
    ## computed pCCA$rank was needed before, but zero it here
    if (!is.null(sol$pCCA) && sol$pCCA$tot.chi == 0)
        pCCA$rank <- 0
    sol <- c(list(call = call), sol)
    sol$method <- "cca"
    sol$inertia <- "mean squared contingency coefficient"
    class(sol) <- "cca"
    sol
}
