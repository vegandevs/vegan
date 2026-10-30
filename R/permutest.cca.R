permutest <- function(x, ...)
    UseMethod("permutest")

permutest.default <- function(x, ...)
    stop("no default permutation test defined")

`permutest.cca` <-
    function (x, permutations = how(nperm=99),
              model = c("reduced", "direct", "full"), by = NULL, first = FALSE,
              strata = NULL, parallel = getOption("mc.cores"), ...,
              permustat = c("pseudo-F", "modified-F"), groups = NULL)
{
    permustat <- match.arg(permustat)
    if (permustat == "modified-F") {
        if (inherits(x, "capscale"))
            stop("'permustat = \"modified-F\"' is not available for capscale")
        if (!inherits(x, "rda"))
            stop("'permustat = \"modified-F\"' is not available for cca; ",
                 "use adonis2() or dbrda() with explicit dissimilarities")
    }
    ## do something sensible with insensible input (no constraints)
    if (is.null(x$CCA) || x$CCA$rank == 0) {
        if (permustat == "modified-F")
            stop("the model must have a saturated group constraint space ",
                 "for 'permustat = \"modified-F\"'")
        sol <- list(call = match.call(), testcall = x$call, model = NA,
                    F.0 = NA, F.perm = NA, chi = c(0, max(x$CA$tot.chi,0)),
                    num = 0, den = x$CA$tot.chi,
                    df = c(0, nobs(x) - max(x$pCCA$QR$rank,0) - 1),
                    nperm = 0, method = x$method, first = FALSE,
                    termlabels = "<none>", Random.seed = NA)
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
    w <- attr(ordiYbar(x, "initial"), "RW")
    isCCA <- !is.null(w)
    isPartial <- !is.null(x$pCCA) && x$pCCA$rank > 0   # handle conditions
    isDB <- inherits(x, c("dbrda")) # only dbrda is distance-based
    if (permustat == "modified-F") {
        if (!is.null(by))
            stop("'by' must be NULL with 'permustat = \"modified-F\"'")
        if (first)
            stop("'first' must be FALSE with 'permustat = \"modified-F\"'")
        if (model == "full")
            stop("'model = \"full\"' is not available with modified-F")
        group.info <- modifiedFGroups(x, groups)
        Q <- x$CCA$QR
        E <- ordiYbar(x, "initial")
        N <- nrow(E)
        permutations <- getPermuteMatrix(permutations, N, strata = strata)
        nperm <- nrow(permutations)
        observed <- modifiedFStatistic(E, Q, group.info, isDB)

        getF2lapply <- function(i, permutations, E, Q, group.info, isDB) {
            modifiedFPermutations(permutations[i, , drop = FALSE], E, Q,
                                  group.info, isDB)
        }
        if (is.null(parallel))
            parallel <- 1
        hasClus <- inherits(parallel, "cluster")
        if (hasClus || parallel > 1) {
            if (.Platform$OS.type == "unix" && !hasClus) {
                tmp <- do.call(rbind,
                    mclapply(splitIndices(nperm, parallel), getF2lapply,
                             mc.cores = parallel, permutations = permutations,
                             E = E, Q = Q, group.info = group.info,
                             isDB = isDB))
            } else {
                if (!hasClus) {
                    parallel <- makeCluster(parallel)
                    on.exit(stopCluster(parallel))
                }
                tmp <- do.call(rbind,
                    parLapply(parallel, splitIndices(nperm, length(parallel)),
                              getF2lapply, permutations = permutations,
                              E = E, Q = Q, group.info = group.info,
                              isDB = isDB))
            }
        } else {
            tmp <- modifiedFPermutations(permutations, E, Q, group.info, isDB)
        }

        Call <- match.call()
        Call[[1L]] <- as.name("permutest")
        res.df <- nobs(x) - Q$rank - 1
        sol <- list(call = Call, testcall = x$call, model = model,
                    permustat = permustat,
                    F.0 = structure(observed$statistic, names = "Model"),
                    F.perm = matrix(tmp[, "F"], ncol = 1L),
                    chi = c(Model = x$CCA$tot.chi,
                            Residual = x$CA$tot.chi),
                    num = tmp[, "num"], den = tmp[, "den"],
                    den.0 = observed$denominator, den.df = observed$den.df,
                    res.df = res.df,
                    df = c(group.info$rank, observed$den.df),
                    groups = group.info$groups,
                    group.size = group.info$size,
                    group.dispersion = observed$dispersion,
                    dispersion.coefficients = observed$coefficients,
                    nperm = nperm, method = x$method, first = FALSE,
                    termlabels = "Model")
        sol$Random.seed <- attr(permutations, "seed")
        sol$control <- attr(permutations, "control")
        if (!missing(strata)) {
            sol$strata <- deparse1(substitute(strata))
            sol$stratum.values <- strata
        }
        class(sol) <- "permutest.cca"
        return(sol)
    }
    ## C function to get the statistics in one loop
    getF <- function(indx, E, Q, QZ, effects, w, first, isPartial, isCCA,
                     isDB, q, r)
    {
        # q is the rank(s) of the effect(s) - 1
        # r is the rank of the residual term - 1
        if (!is.matrix(indx))
            indx <- matrix(indx, nrow = 1)
        out <- .Call(do_getF, indx, E, Q, QZ, effects, w, first, isPartial,
                     isCCA, isDB, PACKAGE = "vegan")
        p <- length(effects)
        if (!isPartial && !first)
            out[, p + 1] <- Chi.tot - rowSums(out[,seq_len(p), drop=FALSE])
        if (p > 1) {
            if (by == "terms")
                out[, seq_len(p)] <- sweep(out[, seq_len(p), drop = FALSE],
                                           2, q, "/")
            out <- cbind(out, sweep(out[,seq_len(p), drop=FALSE], 1,
                                    out[, p + 1] / r, "/"))
        }
        else
            out <- cbind(out, (out[, 1] / q) / (out[, 2] / r))
        out
    }
    ## end getF
    ## wrapper to getF for in parallel lapply
    getFlapply <- function(i, permutations, E, Q, QZ, effects, w, first,
                            isPartial, isCCA, isDB, q, r) {
        getF(permutations[i,], E = E, Q = Q, QZ = QZ, effects = effects, w = w,
             first = first, isPartial = isPartial, isCCA = isCCA, isDB = isDB,
             q = q, r = r)
    }
    ## wrapper to getF for parRapply
    getFparRapply <- function(i, permutations, E, Q, QZ, effects, w, first,
                              isPartial, isCCA, isDB, q, r) {
        getF(permutations[i,], E = E, Q = Q, QZ = QZ, effects = effects, w = w,
             first = first, isPartial = isPartial, isCCA = isCCA, isDB = isDB,
             q = q, r = r)
    }
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
    ## residualizing constraints will remove pCCA terms
    if (isPartial)
        effects <- effects - x$pCCA$QR$rank

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
                           mclapply(splitIndices(nperm, parallel),
                                    getFlapply,
                                    mc.cores = parallel,
                                    permutations = permutations,
                                    E = E, Q = Q, QZ = QZ,
                                    effects = effects, w = w,
                                    first = first, isPartial = isPartial,
                                    isCCA = isCCA, isDB = isDB, q = q, r = r))
        } else {
            ## if hasClus, do not set up and stop a temporary cluster
            if (!hasClus) {
                parallel <- makeCluster(parallel)
                on.exit(stopCluster(parallel))
            }
            tmp <- do.call(rbind,
                           parLapply(parallel,
                                     splitIndices(nperm, length(parallel)),
                                     getFparRapply,
                                     permutations=permutations,
                                     E = E, Q = Q, QZ = QZ, effects = effects, w = w,
                                     first = first, isPartial = isPartial,
                                     isCCA = isCCA, isDB = isDB, q = q, r = r))
        }
    } else {
        tmp <- getF(permutations, E = E, Q = Q, QZ = QZ, effects = effects,
                    w = w, first = first, isPartial = isPartial, isCCA = isCCA,
                    isDB = isDB, q = q, r = r)
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
                method = x$method, inertia = x$inertia, first = first,
                termlabels = termlabs)
    sol$Random.seed <- attr(permutations, "seed")
    sol$control <- attr(permutations, "control")
    if (!missing(strata)) {
        sol$strata <- deparse1(substitute(strata))
        sol$stratum.values <- strata
    }
    class(sol) <- "permutest.cca"
    sol
}
