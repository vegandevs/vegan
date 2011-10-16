"mantel" <-
  function (xdis, ydis, method = "pearson", permutations = 999, 
            strata) 
{
    xdis <- as.dist(xdis)
    ydis <- as.vector(as.dist(ydis))
    tmp <- cor.test(as.vector(xdis), ydis, method = method)
    statistic <- as.numeric(tmp$estimate)
    variant <- tmp$method
    N <- attr(xdis, "Size")
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
    if (permutations) {
        perm <- numeric(permutations)
        ## asdist as an index selects lower diagonal like as.dist,
        ## but is faster since it does not set 'dist' attributes
        xmat <- as.matrix(xdis)
        asdist <- row(xmat) > col(xmat)
        ptest <- function(take, ...) {
            permvec <- (xmat[take, take])[asdist]
            drop(cor(permvec, ydis, method = method))
        }
        perm <- sapply(1:permutations, function(i, ...) ptest(permat[i,], ...) )
        signif <- (sum(perm >= statistic) + 1)/(permutations + 1)
     }
    else {
        signif <- NA
        perm <- NULL
    }
    res <- list(call = match.call(), method = variant, statistic = statistic, 
                signif = signif, perm = perm, permutations = permutations)
    if (!missing(strata)) {
        res$strata <- deparse(substitute(strata))
        res$stratum.values <- strata
    }
    class(res) <- "mantel"
    res
}
