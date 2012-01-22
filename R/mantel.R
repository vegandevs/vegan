`mantel` <-
  function (xdis, ydis, method = "pearson", permutations = 999, 
            strata, parallel = getOption("mc.cores")) 
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
        ## Parallel processing
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
