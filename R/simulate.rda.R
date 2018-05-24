`simulate.rda` <-
    function(object, nsim = 1, seed = NULL, indx = NULL, rank = "full",
             correlated = FALSE, ...)
{
    ## Fail if there is no constrained component (it could be possible
    ## to change the function to handle unconstrained ordination, too,
    ## when rank < "full", but that would require redesign)
    if (is.null(object$CCA))
        stop("function can be used only with constrained ordination")

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
    ## indx can be an output of permute::shuffleSet in which case it
    ## is a nsim x nrow matrix, or it may be a single vector, in which
    ## case it will changed to shuffleSet
    if (!is.null(indx))
        if (is.vector(indx))
            dim(indx) <- c(1, length(indx))
    ## If nsim is missing, take it from indx (if given)
    if (missing(nsim) && !is.null(indx))
        nsim <- nrow(indx)
    ## Check that dims match
    if (!is.null(indx))
        if(nrow(indx) != nsim)
            stop(gettextf("'nsim' (%d) and no. of 'indx' rows (%d) do not match",
                          nsim, nrow(indx)))
    ## collect data to back-transform data to the scale of observations
    sqnr1 <- sqrt(nobs(object) - 1)
    ## the ifs are only needed to cope with pre-2.5-0 vegan: now
    ## we always have Ybar, but earlier we needed to check whether
    ## we had CA or CCA Xbar
    if (!is.null(object$Ybar)) {
        cnt <- attr(object$Ybar, "scaled:center")
        scl <- attr(object$Ybar, "scaled:scale")
    } else { # needed for vegan-2.4 compatibility
        if (is.null(object$CCA))
            tmp <- object$CA$Xbar
        else tmp <- object$CCA$Xbar
        cnt <- attr(tmp, "scaled:center")
        scl <- attr(tmp, "scaled:scale")
    }

    ## Proper simulation: very similar for simulate.lm, but produces
    ## an array of response matrices

    ftd <- predict(object, type = "working", rank = rank)
    ## pRDA: add partial Fit to the constrained
    if (!is.null(object$pCCA))
        ftd <- ftd + ordiYbar(object, "pCCA")
    ## if(is.null(indx)), we have parametric Gaussian simulation and
    ## need to generate sd matrices. The residuals sd is always taken

    ## from the unconstrained (residual) component. If
    ## species are uncorrelated, we need only species sd's, but if
    ## correlated, we also need species covariances.
    CAYbar <- ordiYbar(object, "CA")
    if (!correlated)
        dev <- outer(rep(1, nrow(ftd)), apply(CAYbar, 2, sd))
    else
        dev <- cov(CAYbar)
    ## Generate an array
    ans <- array(0, c(dim(ftd), nsim))
    for (i in seq_len(nsim)) {
        if (!is.null(indx))
            ans[,,i] <- as.matrix(ftd + CAYbar[indx[i,],])
        else if (!correlated)
            ans[,,i] <- as.matrix(ftd + matrix(rnorm(length(ftd), sd = dev),
                                               nrow = nrow(ftd)))
        else {
            ans[,,i] <- t(apply(ftd, 1,
                                function(x) mvrnorm(1, mu = x, Sigma = dev)))
        }
        ## ans to the scale of observations
        ans[,,i] <- ans[,,i] * sqnr1
        if (!is.null(scl))
            ans[,,i] <- sweep(ans[,,i], 2, scl, "*")
        ans[,,i] <- sweep(ans[,,i], 2, cnt, "+")
    }
    ## set RNG attributes
    if (is.null(indx))
        attr(ans, "seed") <- RNGstate
    else
        attr(ans, "seed") <- "index"
    ## set commsim attributes if nsim > 1, else return a 2-dim matrix
    if (nsim == 1) {
        ans <- ans[,,1]
        attributes(ans) <- attributes(ftd)
    } else {
        dimnames(ans) <- list(rownames(ftd), colnames(ftd),
                              paste("sim", seq_len(nsim), sep = "_"))
        orig <- (ftd + CAYbar) * sqnr1
        if (!is.null(scl))
            orig <- sweep(orig, 2, scl, "*")
        orig <- sweep(orig, 2, cnt, "+")
        attr(ans, "data") <- round(orig, 12)
        attr(ans, "method") <- paste("simulate", ifelse(is.null(indx),
                                                        "parametric", "index"))
        attr(ans, "binary") <- FALSE
        attr(ans, "isSeq") <- FALSE
        attr(ans, "mode") <- "double"
        attr(ans, "start") <- 1L
        attr(ans, "end") <- as.integer(nsim)
        attr(ans, "thin") <- 1L
        class(ans) <- c("simulate.rda", "simmat", "array")
    }
    ans
}

### simulate. cca was cloned from simulate.rda.  Works with internal
### Chi-square standardized form, and at the end back-standardizes
### with row and column totals and matrix grand totals. This does not
### still guarantee that all marginal totals are positive.

`simulate.cca` <-
    function(object, nsim = 1, seed = NULL, indx = NULL, rank = "full",
             correlated = FALSE, ...)
{
    ## Fail if no CCA
    if (is.null(object$CCA))
        stop("function can be used only with constrained ordination")
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
    ## Preparations like in simulate.rda()
    if (!is.null(indx))
        if (is.vector(indx))
            dim(indx) <- c(1, length(indx))
    if (missing(nsim) && !is.null(indx))
        nsim <- nrow(indx)
    if (!is.null(indx))
        if(nrow(indx) != nsim)
            stop(gettextf("'nsim' (%d) and no. of 'indx' rows (%d) do not match",
                          nsim, nrow(indx)))
    ## Need sqrt of rowsums for weighting
    sq.r <- sqrt(object$rowsum)
    ## Fitted value
    ftd <- predict(object, type = "working", rank = rank)
    ## pCCA: add partial Fit to the constrained
    if (!is.null(object$pCCA))
        ftd <- ftd + ordiYbar(object, "pCCA")
    ## Residual Xbar need weighting and back-weighting
    Xbar <- sweep(ordiYbar(object, "CA"), 1, sq.r, "*")
    ## Simulation
    if (correlated)
        dev <- cov(Xbar)
    else
        dev <- outer(rep(1, nrow(ftd)), apply(Xbar, 2, sd))
    ans <- array(0, c(dim(ftd), nsim))
    for (i in seq_len(nsim)) {
        if (is.null(indx)) {
            if (correlated)
                tmp <- mvrnorm(nrow(ftd), numeric(ncol(ftd)), Sigma = dev)
            else
                tmp <- matrix(rnorm(length(ftd), sd = dev),
                          nrow = nrow(ftd))
            ans[,,i] <- as.matrix(ftd + sweep(tmp, 1, sq.r, "/"))
        }
        else
            ans[,,i] <- as.matrix(ftd + sweep(Xbar[indx[i,],], 1, sq.r, "/"))
    }
    ## From internal form to the original form with fixed marginal totals
    rc <- object$rowsum %o% object$colsum
    for (i in seq_len(nsim))
        ans[,,i] <- (ans[,,i] * sqrt(rc) + rc) * object$grand.total
    ## RNG attributes
    if (is.null(indx))
        attr(ans, "seed") <- RNGstate
    else
        attr(ans, "seed") <- "index"
    ## set commsim attributes if nsim > 1, else return a 2-dim matrix
    if (nsim == 1) {
        ans <- ans[,,1]
        attributes(ans) <- attributes(ftd)
    } else {
        dimnames(ans) <- list(rownames(ftd), colnames(ftd),
                              paste("sim", seq_len(nsim), sep = "_"))
        obsdata <- ordiYbar(object, "initial")
        obsdata <- (obsdata * sqrt(rc) + rc) * object$grand.total
        attr(ans, "data") <- round(obsdata, 12)
        attr(ans, "method") <- paste("simulate", ifelse(is.null(indx),
                                                        "parametric", "index"))
        attr(ans, "binary") <- FALSE
        attr(ans, "isSeq") <- FALSE
        attr(ans, "mode") <- "double"
        attr(ans, "start") <- 1L
        attr(ans, "end") <- as.integer(nsim)
        attr(ans, "thin") <- 1L
        class(ans) <- c("simulate.cca", "simmat", "array")
    }
    ans
}


### capscale method: copies simulate.rda as much as possible. Function
### works with the internal metric scaling mapping of fit and error,
### but returns Euclidean distances adjusted to the original scaling
### of input dissimilarities. Only the real components are used, and
### capscale() of simulated dissimilarities have no Imaginary
### component.

`simulate.capscale` <-
    function(object, nsim = 1, seed = NULL, indx = NULL, rank = "full",
             correlated = FALSE, ...)
{
    ## Fail if no CCA component
    if (is.null(object$CCA))
        stop("function can be used only with constrained ordination")
    if (is.null(indx) && correlated)
        warning("argument 'correlated' does not work and will be ignored")
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
        ftd <- ordiYbar(object, "CCA")
        ## redo analysis when rank < full
        if (rank < object$CCA$rank) {
            x <- svd(ftd, nu = rank, nv = rank)
            ftd <- x$u %*% diag(x$d[1:rank], nrow=rank) %*% t(x$v)
        }
    } else {
        ftd <- 0
    }
    ## add partial Fit to the constrained
    if (!is.null(object$pCCA))
        ftd <- ftd + ordiYbar(object, "pCCA")
    if (is.null(indx))
        ans <- as.data.frame(ftd + matrix(rnorm(length(ftd),
               sd = outer(rep(1,nrow(ftd)), apply(ordiYbar(object, "CA"), 2, sd))),
               nrow = nrow(ftd)))
    else
        ans <- ftd + ordiYbar(object, "CA")[indx,]
    ## return Euclidean distances
    ans <- ans * object$adjust
    ans <- dist(ans)
    ## remove adjustment done in capscale and put dissimilarities to
    ## (approximately) original scale
    if (is.null(indx))
        attr(ans, "seed") <- RNGstate
    else
        attr(ans, "seed") <- indx
    ans
}

### simulate.dbrda cannot be done along similar lines as
### simulate.capscale, because low-rank approximation needs column
### scores v and cannot be found only from row scores u that are the
### only ones we have in dbrda(). Residuals also need exra thinking,
### and therefore we just disable simulate.dbrda()

`simulate.dbrda` <-
    function(object, nsim = 1, seed = NULL, ...)
{
    .NotYetImplemented()
}
