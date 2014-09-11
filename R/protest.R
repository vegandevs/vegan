`protest` <-
    function (X, Y, scores = "sites", permutations = how(nperm = 999),
              ...)
{
    X <- scores(X, display = scores, ...)
    Y <- scores(Y, display = scores, ...)
    ## Centre and normalize X & Y here so that the permutations will
    ## be faster
    X <- scale(X, scale = FALSE)
    Y <- scale(Y, scale = FALSE)
    X <- X/sqrt(sum(X^2))
    Y <- Y/sqrt(sum(Y^2))
    ## Transformed X and Y will yield symmetric procrustes() and we
    ## need not specify that in the call (but we set it symmetric
    ## after the call).
    sol <- procrustes(X, Y, symmetric = FALSE)
    sol$symmetric <- TRUE
    sol$t0 <- sqrt(1 - sol$ss)
    N <- nrow(X)

    ## Permutations: We only need the goodness of fit statistic from
    ## Procrustes analysis, and therefore we only have the necessary
    ## function here. This avoids a lot of overhead of calling
    ## procrustes() for each permutation. The following gives the
    ## Procrustes r directly.
    procr <- function(X, Y) sum(svd(crossprod(X, Y), nv=0, nu=0)$d)

    permutations <- getPermuteMatrix(permutations, N)
    if (ncol(permutations) != N)
        stop(gettextf("'permutations' have %d columns, but data have %d observations",
                      ncol(permutations), N))
    np <- nrow(permutations)

    perm <- sapply(seq_len(np),
                   function(i, ...) procr(X, Y[permutations[i,],]))

    Pval <- (sum(perm >= sol$t0) + 1)/(np + 1)

    sol$t <- perm
    sol$signif <- Pval
    sol$permutations <- np
    sol$control <- attr(permutations, "control")
    sol$call <- match.call()
    class(sol) <- c("protest", "procrustes")
    sol
}
