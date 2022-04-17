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

    ## Function will give extremely high values when all species occur
    ## only once or S==N, starting from fisherfit(1) which is ca. 1e8,
    ## and it can make sense to have special treatment of S==N. With S
    ## == 0, we force alpha 0 whereas the function would give
    ## fisherfit(0) as 1 (which hardly makes sense).
    if (S > 0) {
        sol <- uniroot(d1fun, c(1,50), extendInt = "upX", S = S, N = N, ...)
        if (S == N)
            warning("all species singletons: alpha arbitrarily high")
    } else {
        sol <- list(root = 0, iter = 0, estim.prec = NA)
    }
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
