simulate.nullmodel <-
function(object, nsim=1, seed = NULL, burnin=0, thin=1, ...)
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
    if (nsim < 1)
        stop("'nsim' must be at least 1")
    m <- object$data
    if (object$commsim$isSeq) {
        ## here is burnin, see update method
        if (burnin > 0)
            object <- update(object, burnin, ...)
        x <- object$state
    } else {
        x <- m
#        if (thin != 1)
#            message("non-sequential model: 'thin' set to 1")
        thin <- 1L
#        if (burnin != 0)
#            message("non-sequential model: 'burnin' set to 0")
        burnin <- 0L
    }
    perm <- object$commsim$fun(x=x,
        n=as.integer(nsim),
        nr=object$nrow,
        nc=object$ncol,
        rs=object$rowSums,
        cs=object$colSums,
        rf=object$rowFreq,
        cf=object$colFreq,
        s=object$totalSum,
        fill=object$fill,
        thin=as.integer(thin), ...)
    if (object$commsim$isSeq) {
        Start <- as.integer(object$iter + 1L)
        End <- as.integer(object$iter + nsim * thin)
        state <- perm[,,nsim]
        storage.mode(state) <- object$commsim$mode
        assign("state", state, envir=object)
        assign("iter", as.integer(End), envir=object)
    } else {
        Start <- 1L
        End <- as.integer(nsim)
    }
    attr(perm, "data") <- m
    attr(perm, "seed") <- RNGstate
    attr(perm, "method") <- object$commsim$method
    attr(perm, "binary") <- object$commsim$binary
    attr(perm, "isSeq") <- object$commsim$isSeq
    attr(perm, "mode") <- object$commsim$mode
    attr(perm, "start") <- Start
    attr(perm, "end") <- End
    attr(perm, "thin") <- as.integer(thin)
    class(perm) <- c("simmat", "array")
    dimnames(perm) <- list(rownames(m), colnames(m),
        paste("sim", seq_len(nsim), sep = "_"))
    perm
}
