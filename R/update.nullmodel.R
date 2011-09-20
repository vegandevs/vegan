update.nullmodel <-
function(object, nsim=1, seed = NULL, ...)
{
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
    m <- object$data
    if (object$commsim$isSeq) {
        perm <- object$commsim$fun(x=object$state,
            n=nsim,
            nr=object$nrow,
            nc=object$ncol,
            rs=object$rowSums,
            cs=object$colSums,
            rf=object$rowFreq,
            cf=object$colFreq,
            s=object$totalSum,
            fill=object$fill,
            thin=1, ...)
        state <- perm[,,nsim]
        storage.mode(state) <- object$commsim$mode
        assign("state", state, envir=object)
        assign("iter", as.integer(object$iter + nsim), envir=object)
    }
    invisible(NULL)
}
