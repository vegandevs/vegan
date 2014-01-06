`anosim` <-
    function (dat, grouping, permutations = 999,
              distance = "bray", strata, parallel = getOption("mc.cores")) 
{
    if (inherits(dat, "dist")) 
        x <- dat
    else if (is.matrix(dat) && nrow(dat) == ncol(dat) && all(dat[lower.tri(dat)] == 
        t(dat)[lower.tri(dat)])) {
        x <- dat
        attr(x, "method") <- "user supplied square matrix"
    }
    else x <- vegdist(dat, method = distance)
    if (any(x < -sqrt(.Machine$double.eps)))
        warning("some dissimilarities are negative -- is this intentional?")
    sol <- c(call = match.call())
    grouping <- as.factor(grouping)
    matched <- function(irow, icol, grouping) {
        grouping[irow] == grouping[icol]
    }
    x.rank <- rank(x)
    N <- attr(x, "Size")
    div <- length(x)/2
    irow <- as.vector(as.dist(row(matrix(nrow = N, ncol = N))))
    icol <- as.vector(as.dist(col(matrix(nrow = N, ncol = N))))
    within <- matched(irow, icol, grouping)
    aver <- tapply(x.rank, within, mean)
    statistic <- -diff(aver)/div
    cl.vec <- rep("Between", length(x))
    take <- as.numeric(irow[within])
    cl.vec[within] <- levels(grouping)[grouping[take]]
    cl.vec <- factor(cl.vec, levels = c("Between", levels(grouping)))
    ptest <- function(take, ...) {
        cl.perm <- grouping[take]
        tmp.within <- matched(irow, icol, cl.perm)
        tmp.ave <- tapply(x.rank, tmp.within, mean)
        -diff(tmp.ave)/div
    }
    if (length(permutations) == 1) {
        if (permutations > 0) {
            arg <- if (missing(strata)) NULL else strata
            permat <- t(replicate(permutations, permuted.index(N, strata = arg)))
        }
    } else {
        permat <- as.matrix(permutations)
        if (ncol(permat) != N)
            stop(gettextf("'permutations' have %d columns, but data have %d rows",
                          ncol(permat), N))
        permutations <- nrow(permat)
    }
    ## Parallel processing
    if (is.null(parallel))
        parallel <- 1
    hasClus <- inherits(parallel, "cluster")
    if ((hasClus || parallel > 1)  && require(parallel)) {
        if(.Platform$OS.type == "unix" && !hasClus) {
            perm <- unlist(mclapply(1:permutations, function(i, ...)
                                    ptest(permat[i,]),
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
        perm <- sapply(1:permutations, function(i) ptest(permat[i,]))
    }
    p.val <- (1 + sum(perm >= statistic))/(1 + permutations)
    sol$signif <- p.val
    sol$perm <- perm
    sol$permutations <- permutations
    sol$statistic <- as.numeric(statistic)
    sol$class.vec <- cl.vec
    sol$dis.rank <- x.rank
    sol$dissimilarity <- attr(x, "method") 
    if (!missing(strata)) {
        sol$strata <- deparse(substitute(strata))
        sol$stratum.values <- strata
    }
    class(sol) <- "anosim"
    sol
}
