`mantel.partial` <-
  function (xdis, ydis, zdis, method = "pearson", permutations = 999,
            strata = NULL, na.rm = FALSE, parallel = getOption("mc.cores"))
{
    EPS <- sqrt(.Machine$double.eps)
    part.cor <- function(rxy, rxz, ryz) {
        (rxy - rxz * ryz)/sqrt(1-rxz*rxz)/sqrt(1-ryz*ryz)
    }
    if ((is.matrix(xdis) || is.data.frame(xdis)) &&
               isSymmetric(unname(as.matrix(xdis))))
        xdis <- as.dist(xdis)
    if (!inherits(xdis, "dist"))
        stop("'xdis' must be a 'dist' object or a symmetric square matrix")
    if ((is.matrix(ydis) || is.data.frame(ydis)) &&
               isSymmetric(unname(as.matrix(ydis))))
        ydis <- as.dist(ydis)
    if (!inherits(ydis, "dist"))
        stop("'ydis' must be a 'dist' object or a symmetric square matrix")
    if ((is.matrix(zdis) || is.data.frame(zdis)) &&
               isSymmetric(unname(as.matrix(zdis))))
        zdis <- as.dist(zdis)
    if (!inherits(zdis, "dist"))
        stop("'zdis' must be a 'dist' object or a symmetric square matrix")

        ydis <- as.dist(ydis)
    xdis <- as.dist(xdis)
    ydis <- as.vector(ydis)
    zdis <- as.vector(zdis)
    ## Handle missing values
    if (na.rm)
        use <- "complete.obs"
    else
        use <- "all.obs"
    rxy <- cor(as.vector(xdis), ydis, method = method, use = use)
    rxz <- cor(as.vector(xdis), zdis, method = method, use = use)
    ryz <- cor(ydis, zdis, method = method, use = use)
    variant <- match.arg(method, eval(formals(cor)$method))
    variant <- switch(variant,
                      pearson = "Pearson's product-moment correlation",
                      kendall = "Kendall's rank correlation tau",
                      spearman = "Spearman's rank correlation rho",
                      variant)
    statistic <- part.cor(rxy, rxz, ryz)
    N <- attr(xdis, "Size")
    permat <- getPermuteMatrix(permutations, N, strata = strata)
    if (ncol(permat) != N)
        stop(gettextf("'permutations' have %d columns, but data have %d observations",
                      ncol(permat), N))
    permutations <- nrow(permat)

    if (permutations) {
        N <- attr(xdis, "Size")
        perm <- rep(0, permutations)
        xmat <- as.matrix(xdis)
        asdist <- row(xmat) > col(xmat)
        ptest <- function(take, ...) {
            permvec <- (xmat[take, take])[asdist]
            rxy <- cor(permvec, ydis, method = method, use = use)
            rxz <- cor(permvec, zdis, method = method, use = use)
            part.cor(rxy, rxz, ryz)
        }
        ## parallel processing
        if (is.null(parallel))
            parallel <- 1
        hasClus <- inherits(parallel, "cluster")
        if (hasClus || parallel > 1) {
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
        signif <- (sum(perm >= statistic - EPS)+1)/(permutations + 1)
    }
    else {
        signif <- NA
        perm <- NULL
    }
    res <- list(call = match.call(), method = variant, statistic = statistic,
                signif = signif, perm = perm, permutations = permutations,
                control = attr(permat, "control"))
    if (!missing(strata)) {
        res$strata <- deparse(substitute(strata))
        res$stratum.values <- strata
    }
    class(res) <- c("mantel.partial", "mantel")
    res
}

