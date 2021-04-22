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

    ## need at least 2 species to estimate Fisher alpha -- set to NA
    ## (or 1?) for zero-species and one-species communities
    if (S > 1)
        sol <- uniroot(d1fun, c(1,50), extendInt = "upX", S = S, N = N, ...)
    else
        sol <- list(root = NA, iter = 0, estim.prec = NA)
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
