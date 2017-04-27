## Fisher alpha is actually based only on the number of species S and
## number of individuals.

`fisherfit` <-
    function(x, ...)
{
    nr <- as.fisher(x)
    S <- sum(nr)
    N <- sum(x)
    ## Solve 'x' (Fisher alpha).
    d1fun <- function(x, S, N) x * log(1 + N/x) - S

    ## 'extendInt' arg was added in R r63162 | maechler | 2013-07-03
    ## 11:47:22 +0300 (Wed, 03 Jul 2013) and released in R 3.1.0
    ## (2014-04-10).

    sol <- uniroot(d1fun, c(1,50), extendInt = "upX", S = S, N = N, ...)

    nuisance <- N/(N + sol$root)
    ## we used nlm() earlier, and the following output is compatible
    out <- list(estimate = sol$root, hessian = NA,
                iterations = sol$iter, df.residual = NA,
                nuisance = nuisance, fisher = nr,
                estim.prec = sol$estim.prec,
                code = 2*is.na(sol$estim.prec) + 1)
    class(out) <- "fisherfit"
    out
}
