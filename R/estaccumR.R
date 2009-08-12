##" Individual based accumulation model. Similar to poolaccum but uses
##estimateR. Inherits from "poolaccum" class and uses its methods.
`estaccumR` <-
    function(x, permutations = 100)
{
    n <- nrow(x)
    N <- seq_len(n)
    S <- chao <- ace <- matrix(0, nrow = n, ncol = permutations)
    for (i in 1:permutations) {
        take <- sample(n)
        tmp <- estimateR(apply(x[take,], 2, cumsum))
        S[,i] <- tmp[1,]
        chao[,i] <- tmp[2,]
        ace[, i] <- tmp[4,]
    }
    means <- cbind(N = N, S = rowMeans(S), Chao = rowMeans(chao),
                   ACE = rowMeans(ace))
    out <- list(S = S, chao = chao, ace = ace, N = N, means = means)
    class(out) <- c("estaccumR", "poolaccum")
    out
}
