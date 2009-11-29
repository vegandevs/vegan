`simulate.rda` <-
    function(object, nsim = 1, seed = NULL, indx = NULL, ...) 
{
    ## First check cases that won't work (yet?)
    if (!is.null(object$pCCA))
        stop("not yet implemented for partial models")
    ## Handle RNG: code directly from stats::simulate.lm
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
        runif(1)
    if (is.null(seed)) 
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    ## Proper simulation: very similar for simulate.lm, but produces a
    ## response matrix.
    if (nsim > 1)
        .NotYetUsed("nsim")
    ftd <- fitted(object)
    if (is.null(indx))
        ans <- as.data.frame(ftd + matrix(rnorm(length(ftd), 
               sd = outer(rep(1,nrow(ftd)), sd(object$CA$Xbar))), 
               nrow = nrow(ftd)))
    else
        ans <- as.data.frame(ftd + object$CA$Xbar[indx,])
    attr(ans, "seed") <- RNGstate
    ans
}

`simulate.cca` <-
    function(object, nsim = 1, seed = NULL, ...)
{
    .NotYetImplemented()
}

`simulate.capscale` <-
    function(object, nsim = 1, seed = NULL, ...)
{
    .NotYetImplemented()
}
