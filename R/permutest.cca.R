permutest <- function(x, ...)
    UseMethod("permutest")

permutest.default <- function(x, ...)
    stop("no default permutation test defined")

`permutest.cca` <-
    function (x, permutations = how(nperm=99),
              model = c("reduced", "direct", "full"), by = NULL, first = FALSE,
              strata = NULL, parallel = getOption("mc.cores") ,  ...)
{
    ## do something sensible with insensible input (no constraints)
    if (is.null(x$CCA)) {
        sol <- list(call = match.call(), testcall = x$call, model = NA,
                    F.0 = NA, F.perm = NA, chi = c(0, x$CA$tot.chi),
                    num = 0, den = x$CA$tot.chi,
                    df = c(0, nrow(x$CA$u) - max(x$pCCA$QR$rank,0) - 1),
                    nperm = 0, method = x$method, first = FALSE,
                    Random.seed = NA)
        class(sol) <- "permutest.cca"
        return(sol)
    }
    ## compatible arguments?
    if (!is.null(by)) {
        if (first)
            stop("'by' cannot be used with option 'first=TRUE'")
        by <- match.arg(by, c("onedf", "terms"))
        if (by == "terms" && is.null(x$terminfo))
            stop("by='terms' needs a model fitted with a formula")
    }
    model <- match.arg(model)
    ## special cases
    isCCA <- !inherits(x, "rda")    # weighting
    isPartial <- !is.null(x$pCCA)   # handle conditions
    isDB <- inherits(x, c("dbrda")) # only dbrda is distance-based
    ## C function to get the statististics in one loop
    getF <- function(indx, ...)
    {
        if (!is.matrix(indx))
            indx <- matrix(indx, nrow=1)
        out <- .Call(do_getF, indx, E, Q, QZ, effects, first, isPartial, isDB)
        p <- length(effects)
        if (!isPartial && !first)
            out[,p+1] <- Chi.tot - rowSums(out[,seq_len(p), drop=FALSE])
        if (p > 1) {
            if (by == "terms")
                out[, seq_len(p)] <- sweep(out[, seq_len(p), drop = FALSE],
                                               2, q, "/")
            out <- cbind(out, sweep(out[,seq_len(p), drop=FALSE], 1,
                                    out[,p+1]/r, "/"))
        }
        else
            out <- cbind(out, (out[,1]/q)/(out[,2]/r))
        out
    }
    ## end getF
    ## QR decomposition
        Q <- x$CCA$QR
    if (isPartial) {
        QZ <- x$pCCA$QR
    } else {
        QZ <- NULL
    }
    ## statistics: overall tests
    if (first) {
        Chi.z <- x$CCA$eig[1]
        q <- 1
    }
    else {
        Chi.z <- x$CCA$tot.chi
        names(Chi.z) <- "Model"
        q <- x$CCA$qrank
    }
    ## effects
    if (!is.null(by)) {
        partXbar <- ordiYbar(x, "partial")
        if (by == "onedf") {
            effects <- seq_len(q)
            termlabs <-
                if (isPartial)
                    colnames(Q$qr)[effects + x$pCCA$QR$rank]
                else
                    colnames(Q$qr)[effects]
        } else {                   # by = "terms"
            ass <- x$terminfo$assign
            ## ass was introduced in vegan_2.5-0
            if (is.null(ass))
                stop("update() old ordination result object")
            pivot <- Q$pivot
            if (isPartial)
                pivot <- pivot[pivot > x$pCCA$QR$rank] - x$pCCA$QR$rank
            ass <- ass[pivot[seq_len(x$CCA$qrank)]]
            effects <- cumsum(rle(ass)$length)
            termlabs <- labels(terms(x$terminfo))
            if (isPartial)
                termlabs <- termlabs[termlabs %in% labels(terms(x))]
            termlabs <-termlabs[unique(ass)]
        }
        q <- diff(c(0, effects)) # d.o.f.
        if (isPartial)
            effects <- effects + x$pCCA$QR$rank
        F.0 <- numeric(length(effects))
        for (k in seq_along(effects)) {
            fv <- qr.fitted(Q, partXbar, k = effects[k])
            F.0[k] <- if (isDB) sum(diag(fv)) else sum(fv^2)
        }
    }
    else {
        effects <- 0
        termlabs <- "Model"
    }
    ## Set up
    Chi.xz <- x$CA$tot.chi
    names(Chi.xz) <- "Residual"
    r <- nobs(x) - Q$rank - 1
    if (model == "full")
        Chi.tot <- Chi.xz
    else Chi.tot <- Chi.z + Chi.xz
    if (is.null(by))
        F.0 <- (Chi.z/q)/(Chi.xz/r)
    else {
        Chi.z <- numeric(length(effects))
        for (k in seq_along(effects)) {
            fv <- qr.fitted(Q, partXbar, k = effects[k])
            Chi.z[k] <- if (isDB) sum(diag(fv)) else sum(fv^2)
        }
        Chi.z <- diff(c(0, F.0))
        F.0 <- Chi.z/q * r/Chi.xz
    }

    ## permutation data
    E <- switch(model,
                "direct" = ordiYbar(x, "initial"),
                "reduced" = ordiYbar(x, "partial"),
                "full" = ordiYbar(x, "CA"))
    ## vegan < 2.5-0 cannot use direct model in partial dbRDA
    if (is.null(E) && isDB && isPartial)
        stop("'direct' model cannot be used in old partial-dbrda: update() result")

    ## Save dimensions
    N <- nrow(E)
    permutations <- getPermuteMatrix(permutations, N, strata = strata)
    nperm <- nrow(permutations)
    ## Parallel processing (similar as in oecosimu)
    if (is.null(parallel))
        parallel <- 1
    hasClus <- inherits(parallel, "cluster")
    if (hasClus || parallel > 1) {
        if(.Platform$OS.type == "unix" && !hasClus) {
            tmp <- do.call(rbind,
                           mclapply(1:nperm,
                                    function(i) getF(permutations[i,]),
                                    mc.cores = parallel))
        } else {
            ## if hasClus, do not set up and stop a temporary cluster
            if (!hasClus) {
                parallel <- makeCluster(parallel)
            }
            tmp <- parRapply(parallel, permutations, function(i) getF(i))
            tmp <- matrix(tmp, ncol=3, byrow=TRUE)
            if (!hasClus)
                stopCluster(parallel)
        }
    } else {
        tmp <- getF(permutations)
    }
    if ((p <- length(effects)) > 1) {
        num <- tmp[,seq_len(p)]
        den <- tmp[,p+1]
        F.perm <- tmp[, seq_len(p) + p + 1]
    } else {
        num <- tmp[,1]
        den <- tmp[,2]
        F.perm <- tmp[,3, drop=FALSE]
    }
    Call <- match.call()
    Call[[1]] <- as.name("permutest")
    sol <- list(call = Call, testcall = x$call, model = model,
                F.0 = F.0, F.perm = F.perm,  chi = c(Chi.z, Chi.xz),
                num = num, den = den, df = c(q, r), nperm = nperm,
                method = x$method, first = first, termlabels = termlabs)
    sol$Random.seed <- attr(permutations, "seed")
    sol$control <- attr(permutations, "control")
    if (!missing(strata)) {
        sol$strata <- deparse(substitute(strata))
        sol$stratum.values <- strata
    }
    class(sol) <- "permutest.cca"
    sol
}
