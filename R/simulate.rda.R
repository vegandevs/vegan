`simulate.rda` <-
    function(object, nsim = 1, seed = NULL, indx = NULL, rank = "full", ...) 
{
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
    ftd <- predict(object, type = "response", rank = rank)
    ## pRDA: add partial Fit to the constrained
    if (!is.null(object$pCCA))
        ftd <- ftd + object$pCCA$Fit
    if (is.null(indx))
        ans <- as.data.frame(ftd + matrix(rnorm(length(ftd), 
               sd = outer(rep(1,nrow(ftd)), sd(object$CA$Xbar))), 
               nrow = nrow(ftd)))
    else
        ans <- as.data.frame(ftd + object$CA$Xbar[indx,])
    if (is.null(indx))
        attr(ans, "seed") <- RNGstate
    else
        attr(ans, "seed") <- indx
    ans
}

### simulate. cca was cloned from simulate.rda.  Works with internal
### Chi-square standardized form, and at the end back-standardizes
### with row and column totals and matrix grand totals. This does not
### still guarantee that all marginal totals are positive.

`simulate.cca` <-
    function(object, nsim = 1, seed = NULL, indx = NULL, rank = "full", ...)
{
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
    ## Need sqrt of rowsums for weighting
    sq.r <- sqrt(object$rowsum)
    ## Fitted value
    ftd <- predict(object, type = "working", rank = rank)
    ## pCCA: add partial Fit to the constrained
    if (!is.null(object$pCCA))
        ftd <- ftd + object$pCCA$Fit
    ## Residual Xbar need weighting and back-weighting
    Xbar <- sweep(object$CA$Xbar, 1, sq.r, "*")
    if (is.null(indx)) {
        ans <- matrix(rnorm(length(ftd), 
               sd = outer(rep(1,nrow(ftd)), sd(Xbar))), 
               nrow = nrow(ftd))
        ans <- as.data.frame(ftd + sweep(ans, 1, sq.r, "/"))
    }
    else 
        ans <- as.data.frame(ftd + sweep(Xbar[indx,], 1, sq.r, "/"))
    ## From internal form to the original form with fixed marginal totals
    rc <- object$rowsum %o% object$colsum
    ans <- (ans * sqrt(rc) + rc) * object$grand.total
    if (is.null(indx))
        attr(ans, "seed") <- RNGstate
    else
        attr(ans, "seed") <- indx
    ans
}


### capscale method: copies simulate.rda as much as possible. Function
### works with the internal metric scaling mapping of fit and error,
### but returns Euclidean distances adjusted to the original scaling
### of input dissimilarities. Only the real components are used, and
### capscale() of simulated dissimilarities have no Imaginary
### component.

`simulate.capscale` <-
    function(object, nsim = 1, seed = NULL, indx = NULL, rank = "full", ...) 
{
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
    if (nsim > 1)
        .NotYetUsed("nsim")
    ## predict.capscale cannot be used because it returns either
    ## dissimilarities ("response") or scores with the rank of the
    ## constrained solution, and we need rank of the data (not of
    ## constraints).
    if (rank > 0) {
        ftd <- qr.fitted(object$CCA$QR, object$CCA$Xbar)
        ## redo analysis when rank < full
        if (rank < object$CCA$rank) {
            x <- svd(ftd, nu = rank, nv = rank)
            ftd <- x$u %*% diag(x$d[1:rank], nrow=rank) %*% t(x$v)
        }
    } else {
        ftd <- matrix(0, nrow=nrow(object$CA$Xbar),
                      ncol = ncol(object$CA$Xbar))
    }
    ## add partial Fit to the constrained
    if (!is.null(object$pCCA))
        ftd <- ftd + object$pCCA$Fit
    if (is.null(indx))
        ans <- as.data.frame(ftd + matrix(rnorm(length(ftd), 
               sd = outer(rep(1,nrow(ftd)), sd(object$CA$Xbar))), 
               nrow = nrow(ftd)))
    else
        ans <- ftd + object$CA$Xbar[indx,]
    ## return Euclidean distances
    ans <- dist(ans)
    ## remove adjustment done in capscale and put dissimilarities to
    ## (approximately) original scale
    ans <- ans/object$adjust
    if (is.null(indx))
        attr(ans, "seed") <- RNGstate
    else
        attr(ans, "seed") <- indx
    ans
}

