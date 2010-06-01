"mantel" <-
  function (xdis, ydis, method = "pearson", permutations = 999, 
            strata) 
{
    xdis <- as.dist(xdis)
    ydis <- as.vector(as.dist(ydis))
    tmp <- cor.test(as.vector(xdis), ydis, method = method)
    statistic <- as.numeric(tmp$estimate)
    variant <- tmp$method
    if (permutations) {
        N <- attributes(xdis)$Size
        perm <- rep(0, permutations)
        ## asdist asn an index selects lower diagonal like as.dist,
        ## but is faster since it does not seet 'dist' attributes
        xmat <- as.matrix(xdis)
        asdist <- row(xmat) > col(xmat)
        for (i in 1:permutations) {
            take <- permuted.index(N, strata)
            permvec <- (xmat[take, take])[asdist]
            perm[i] <- cor(permvec, ydis, method = method)
        }
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
