`protest` <-
    function (X, Y, scores = "sites", permutations = 999, strata, ...)
{
    X <- scores(X, display = scores, ...)
    Y <- scores(Y, display = scores, ...)
    sol <- procrustes(X, Y, symmetric = TRUE)
    sol$t0 <- sqrt(1 - sol$ss)
    N <- nrow(X)
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
                   function(i, ...) procrustes(X, Y[permat[i,],],
                                               symmetric = TRUE)$ss)
    perm <- sqrt(1 - perm)
    perm <- c(sol$t0, perm)
    Pval <- sum(perm >= sol$t0)/(permutations + 1)
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
