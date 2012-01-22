`mantel.partial` <-
  function (xdis, ydis, zdis, method = "pearson", permutations = 999, 
            strata, parallel = getOption("mc.cores")) 
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
    N <- attr(xdis, "Size")
    if (length(permutations) == 1) {
        if (permutations > 0) {
            arg <- if(missing(strata)) NULL else strata
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
        N <- attr(xdis, "Size")
        perm <- rep(0, permutations)
        xmat <- as.matrix(xdis)
        asdist <- row(xmat) > col(xmat)
        ptest <- function(take, ...) {
            permvec <- (xmat[take, take])[asdist]
            rxy <- cor(permvec, ydis, method = method)
            rxz <- cor(permvec, zdis, method = method)
            part.cor(rxy, rxz, ryz)
        }
        ## parallel processing
        if (is.null(parallel) && getRversion() >= "2.15.0")
            parallel <- get("default", envir = parallel:::.reg)
        if (is.null(parallel) || getRversion() < "2.14.0")
            parallel <- 1
        hasClus <- inherits(parallel, "cluster")
        if ((hasClus || parallel > 1)  && require(parallel)) {
            if(.Platform$OS.type == "unix" && !hasClus) {
                perm <- do.call(rbind,
                               mclapply(1:permutations,
                                        function(i, ...) ptest(permat[i,],...),
                                        mc.cores = parallel))
            } else {
                if (!hasClus) {
                    parallel <- makeCluster(parallel)
                    clusterEvalQ(parallel, library(vegan))
                }
                perm <- parRapply(parallel, permat, ptest)
                if (!hasClus)
                    stopCluster(parallel)
            }
        } else {
            perm <- sapply(1:permutations, function(i, ...) ptest(permat[i,], ...))
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

