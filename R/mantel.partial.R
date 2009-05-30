"mantel.partial" <-
  function (xdis, ydis, zdis, method = "pearson", permutations = 999, 
            strata) 
{
    part.cor <- function(rxy, rxz, ryz) {
        (rxy - rxz * ryz)/sqrt(1-rxz*rxz)/sqrt(1-ryz*ryz)
    }
    xdis <- as.dist(xdis)
    ydis <- as.vector(as.dist(ydis))
    zdis <- as.vector(as.dist(zdis))
    rxy <- cor.test(as.vector(xdis), ydis, method = method)
    rxz <- cor(as.vector(xdis), zdis, method = method)
    ryz <- cor(ydis, zdis, method = method)
    variant <- rxy$method
    rxy <- rxy$estimate
    statistic <- part.cor(rxy, rxz, ryz)
    if (permutations) {
        N <- attributes(xdis)$Size
        perm <- rep(0, permutations)
        for (i in 1:permutations) {
            take <- permuted.index(N, strata)
            permvec <- as.vector(as.dist(as.matrix(xdis)[take, 
                                                         take]))
            rxy <- cor(permvec, ydis, method = method)
            rxz <- cor(permvec, zdis, method = method)
            perm[i] <- part.cor(rxy, rxz, ryz)
        }
        signif <- (sum(perm >= statistic)+1)/(permutations + 1)
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
    class(res) <- c("mantel.partial", "mantel")
    res
}

