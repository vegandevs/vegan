`mantel` <-
  function (xdis, ydis, method = "pearson", permutations = 999, 
            strata = NULL, na.rm = FALSE, parallel = getOption("mc.cores")) 
{
    xdis <- as.dist(xdis)
    ydis <- as.vector(as.dist(ydis))
    ## Handle missing values
    if (na.rm)
        use <- "complete.obs"
    else
        use <- "all.obs"
    statistic <- cor(as.vector(xdis), ydis, method = method, use = use)
    variant <- match.arg(method, eval(formals(cor)$method))
    variant <- switch(variant,
                      pearson = "Pearson's product-moment correlation",
                      kendall = "Kendall's rank correlation tau",
                      spearman = "Spearman's rank correlation rho",
                      variant)
    N <- attr(xdis, "Size")
    permat <- getPermuteMatrix(permutations, N, strata = strata)
    if (ncol(permat) != N)
        stop(gettextf("'permutations' have %d columns, but data have %d observations",
                      ncol(permat), N))
    permutations <- nrow(permat)

    if (permutations) {
        perm <- numeric(permutations)
        ## asdist as an index selects lower diagonal like as.dist,
        ## but is faster since it does not set 'dist' attributes
        xmat <- as.matrix(xdis)
        asdist <- row(xmat) > col(xmat)
        ptest <- function(take, ...) {
            permvec <- (xmat[take, take])[asdist]
            drop(cor(permvec, ydis, method = method, use = use))
        }
        ## Parallel processing
        if (is.null(parallel))
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
                signif = signif, perm = perm, permutations = permutations,
                control = attr(permat, "control"))
    class(res) <- "mantel"
    res
}
