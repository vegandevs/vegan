`protest` <-
    function (X, Y, scores = "sites", permutations = 999, strata, ...)
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
    
    if (length(permutations) == 1) {
        if (permutations > 0) {
            arg <- if (missing(strata)) NULL else strata
            permat <- t(replicate(permutations,
                                  permuted.index(N, strata = arg)))
        }
    } else {
        permat <- as.matrix(permutations)
        if (ncol(permat) != N)
            stop(gettextf("'permutations' have %d columns, but data have %d observations",
                          ncol(permat), N))
        permutations <- nrow(permutations)
    }
    perm <- sapply(1:permutations,
                   function(i, ...) procr(X, Y[permat[i,],]))
    Pval <- (sum(perm >= sol$t0) + 1)/(permutations + 1)
    if (!missing(strata)) {
        strata <- deparse(substitute(strata))
        s.val <- strata
    }
    else {
        strata <- NULL
        s.val <- NULL
    }
    sol$t <- perm
    sol$signif <- Pval
    sol$permutations <- permutations
    sol$strata <- strata
    sol$stratum.values <- s.val
    sol$call <- match.call()
    class(sol) <- c("protest", "procrustes")
    sol
}
