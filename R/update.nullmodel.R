update.nullmodel <-
function(object, nsim=1, seed = NULL, ...)
{
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1)
    if (!is.null(seed)) {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    if (object$commsim$isSeq) {
        perm <- object$commsim$fun(x=object$state,
            n=1L,
            nr=object$nrow,
            nc=object$ncol,
            rs=object$rowSums,
            cs=object$colSums,
            rf=object$rowFreq,
            cf=object$colFreq,
            s=object$totalSum,
            fill=object$fill,
            thin=as.integer(nsim), ...)
        state <- perm[,,1L]
        storage.mode(state) <- object$commsim$mode
        iter <- as.integer(object$iter + nsim)
#        assign("state", state, envir=object)
#        assign("iter", iter, envir=object)
#        attr(state, "iter") <- iter
        out <- nullmodel(state, object$commsim)
        out$iter <- iter
        out$data <- object$data
    } else {
#        state <- NULL
        out <- object
    }
#    invisible(state)
    out
}
