`poolaccum` <-
    function(x, permutations = 100, minsize = 3)
{
    x <- as.matrix(x)
    n <- nrow(x)
    m <- ncol(x)
    N <- seq_len(n)
    S <- chao <- boot <- jack1 <- jack2 <-
        matrix(0, nrow=n, ncol=permutations)
    ## specpool() is slow, but the vectorized versions below are
    ## pretty fast
    for (i in 1:permutations) {
        take <- sample.int(n, n)
        tmp <- apply(x[take,] > 0, 2, cumsum)
        S[,i] <- rowSums(tmp > 0)
        ## All-zero species are taken as *known* to be missing in
        ## subsamples, and in the following we subtract them (as
        ## 2*S-m) from the bootstrap samples to give a more unbiased
        ## estimate.
        boot[,i] <- 2*S[,i] - m + rowSums(exp(sweep(log1p(-sweep(tmp, 1, N, "/")), 1, N, "*") ))
        a1 <- rowSums(tmp == 1)
        a2 <- rowSums(tmp == 2)
        chao[, i] <- S[,i] + ifelse(a2 > 0, a1*a1/2/a2, 0)
        jack1[,i] <- S[,i] + a1 * (N-1)/N
        jack2[,i] <- S[,i] + a1*(2*N-3)/N - a2*(N-2)^2/N/(N-1)
    }
    means <- cbind(`N` = N, `S` = rowMeans(S), `Chao` =  rowMeans(chao),
                   `Jackknife 1` = rowMeans(jack1),
                   `Jackknife 2` = rowMeans(jack2),
                   `Bootstrap` = rowMeans(boot))
    take <- N >= minsize
    out <- list(S = S[take,], chao = chao[take,], jack1 = jack1[take,],
                jack2 = jack2[take,], boot = boot[take,], N = N[take],
                means = means[take,])
    class(out) <- "poolaccum"
    out
}
