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
    ## We may need to bracket the interval
    hi <- 50
    lo <- 1
    tries <- 0
    repeat {
        sol <- try(uniroot(d1fun, c(lo, hi), S = S, N = N, ...), silent = TRUE)
        if (inherits(sol, "try-error")) {
            if(d1fun(hi, S, N) < 0)
                hi <- 2*hi
            if(d1fun(lo, S, N) > 0)
                lo <- lo/2
            tries <- tries + 1
        }
        else break
        ## alpha can tend to +Inf: set root = NA etc.
        if (tries > 200) {
            sol <- list(root = NA, f.root = NA, iter = NA, init.it = NA,
                        estim.prec = NA)
            break
        }
    }
    ## 'extendInt' arg was added in R r63162 | maechler | 2013-07-03
    ## 11:47:22 +0300 (Wed, 03 Jul 2013). Latest release is R 3.0.2 of
    ## 2013-09-25, but it still does not have the argument.  In the
    ## future we may switch to the following:

    ##sol <- uniroot(d1fun, c(1,50), extendInt = "yes", S = S, N = N, ...)
    
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
