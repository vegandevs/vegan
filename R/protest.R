"protest" <-
    function (X, Y, scores = "sites", permutations = 1000, strata, ...) 
{
    X <- scores(X, display = scores, ...)
    Y <- scores(Y, display = scores, ...)
    sol <- procrustes(X, Y, symmetric = TRUE)
    sol$t0 <- sqrt(1 - sol$ss)
    N <- nrow(X)
    perm <- rep(0, permutations)
    for (i in 1:permutations) {
        take <- permuted.index(N, strata)
        tmp <- procrustes(X, Y[take, ], symmetric = TRUE)$ss
        perm[i] <- sqrt(1 - tmp)
    }
    Pval <- sum(perm >= sol$t0)/permutations
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
