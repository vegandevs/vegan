`cca.default` <-
    function (X, Y = NULL, Z = NULL, ...)
{
    ## Protect against grave misuse: some people have used
    ## dissimilarities instead of data
    if (inherits(X, "dist") || NCOL(X) == NROW(X) &&
        isTRUE(all.equal(X, t(X))))
        stop("function cannot be used with (dis)similarities")
    X <- as.matrix(X)
    if (!is.null(Y)) {
        if (is.data.frame(Y) || is.factor(Y)) { # save Y for centroids
            mframe <- as.data.frame(Y) # can be a single factor
            Y <- model.matrix(~ ., as.data.frame(Y))[,-1,drop=FALSE]
        }
        Y <- as.matrix(Y)
    }
    if (!is.null(Z)) {
        if (is.data.frame(Z) || is.factor(Z))
            Z <- model.matrix(~ ., as.data.frame(Z))[,-1,drop=FALSE]
        Z <- as.matrix(Z)
    }

    if (any(rowSums(X) <= 0))
        stop("all row sums must be >0 in the community data matrix")
    if (any(tmp <- colSums(X) <= 0)) {
        exclude.spec <- seq(along=tmp)[tmp]
        names(exclude.spec) <- colnames(X)[tmp]
        class(exclude.spec) <- "exclude"
        X <- X[, !tmp, drop = FALSE]
    }
    sol <- ordConstrained(X, Y, Z, method = "cca")
    ## mframe exists only if function was called as cca(X, mframe)
    if (exists("mframe"))
        sol$CCA$centroids <- getCentroids(sol, mframe)

    if (exists("exclude.spec")) {
        if (!is.null(sol$CCA$v))
            attr(sol$CCA$v, "na.action") <- exclude.spec
        if (!is.null(sol$CA$v))
            attr(sol$CA$v, "na.action") <- exclude.spec
    }
    call <- match.call()
    call[[1]] <- as.name("cca")
    sol <- c(list(call = call,
                  inertia =  "scaled Chi-square"),
             sol)
    class(sol) <- "cca"
    sol
}
