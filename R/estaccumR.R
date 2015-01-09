##" Individual based accumulation model. Similar to poolaccum but uses
##estimateR. Inherits from "poolaccum" class and uses its methods.
`estaccumR` <-
    function(x, permutations = 100, parallel = getOption("mc.cores"))
{
    n <- nrow(x)
    N <- seq_len(n)
    estFun <- function(idx) {
        estimateR(apply(x[idx,], 2, cumsum))[c(1,2,4),]
    }
    permat <- getPermuteMatrix(permutations, n)
    nperm <- nrow(permat)
    ## parallel processing
    if (is.null(parallel))
        parallel <- 1
    hasClus <- inherits(parallel, "cluster")
    if (hasClus || parallel > 1) {
        if(.Platform$OS.type == "unix" && !hasClus) {
            tmp <- mclapply(1:nperm, function(i)
                            estFun(permat[i,]),
                            mc.cores = parallel)
        } else {
            if (!hasClus) {
                parallel <- makeCluster(parallel)
            }
            tmp <- parLapply(parallel, 1:nperm, function(i) estFun(permat[i,]))
            if (!hasClus)
                stopCluster(parallel)
        }
    } else {
        tmp <- lapply(1:permutations, function(i) estFun(permat[i,]))
    }

    S <- sapply(tmp, function(x) x[1,])
    chao <- sapply(tmp, function(x) x[2,])
    ace <- sapply(tmp, function(x) x[3,])
    means <- cbind(N = N, S = rowMeans(S), Chao = rowMeans(chao),
                   ACE = rowMeans(ace))
    out <- list(S = S, chao = chao, ace = ace, N = N, means = means)
    attr(out, "control") <- attr(permat, "control")
    class(out) <- c("estaccumR", "poolaccum")
    out
}
