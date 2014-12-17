##" Individual based accumulation model. Similar to poolaccum but uses
##estimateR. Inherits from "poolaccum" class and uses its methods.
`estaccumR` <-
    function(x, permutations = 100)
{
    n <- nrow(x)
    N <- seq_len(n)
    estFun <- function(idx) {
        estimateR(apply(x[idx,], 2, cumsum))[c(1,2,4),]
    }
    permat <- getPermuteMatrix(permutations, n)
    tmp <- lapply(1:nperm, function(i) estFun(permat[i,]))
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
