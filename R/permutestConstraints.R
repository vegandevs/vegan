## Subset an unevaluated design to the model's retained rows before generation.
constraintPermutationDesign <- function(x, permutations) {
    if (!inherits(permutations, "how") || is.null(getPlots(permutations)$variables))
        return(permutations)
    plots <- getPlots(permutations)
    data <- plots$variables
    target <- rownames(ordiYbar(x, "initial"))
    if (identical(rownames(data), target)) return(permutations)
    take <- match(target, rownames(data))
    if (anyNA(take)) {
        if (nrow(data) == length(target)) return(permutations)
        stop("plot variables cannot be aligned with the fitted observations")
    }
    plots <- update(plots, strata = getStrata(plots)[take],
                    variables = data[take, , drop = FALSE])
    blocks <- getBlocks(permutations)
    if (!is.null(blocks)) blocks <- blocks[take]
    update(permutations, plots = plots, blocks = blocks, all.perms = NULL)
}

constraintStatistic <- function(response, X, Z) {
    fit <- ordConstrained(response, X, Z, method = "pass")
    q <- if (is.null(fit$CCA)) 0L else fit$CCA$qrank
    r <- nrow(response) - max(fit$CCA$QR$rank, fit$pCCA$QR$rank, 0) - 1L
    num <- if (is.null(fit$CCA)) 0 else fit$CCA$tot.chi
    den <- fit$CA$tot.chi
    if (r < 1L) stop("a constraint permutation leaves no residual degrees of freedom")
    F <- if (q > 0 && r > 0) (num / q) / (den / r) else 0
    c(F = F, num = num, den = den, q = q, r = r)
}

permutestConstraints <- function(x, permutations, strata, parallel, permustat, groups) {
    if (permustat != "pseudo-F")
        stop('model = "constraints" currently supports permustat = "pseudo-F"')
    if (is.null(x$CCA) || x$CCA$qrank < 1L)
        stop("constraint permutations require non-empty constraints")
    permutations <- constraintPermutationDesign(x, permutations)
    response <- ordiYbar(x, "initial")
    N <- nrow(response)
    permutations <- getPermuteMatrix(permutations, N, strata = strata, model = "constraints")
    typed <- inherits(permutations, "plotPermutations")
    indices <- if (typed) permute::permutationIndices(permutations) else as.matrix(permutations)
    if (ncol(indices) != N || anyNA(indices) || any(indices < 1L | indices > N))
        stop("constraint permutation indices must refer to fitted data rows")
    if (!typed && any(apply(indices, 1L, function(i)
        !identical(sort(as.integer(i)), seq_len(N)))))
        stop("ordinary constraint permutations must contain each row index once")
    ## The joint QR stores the centred conditions before the constraints,
    ## even when the conditions explain zero inertia.
    Q <- x$CCA$QR
    QZ <- x$pCCA$QR
    ncond <- length(QZ$pivot)
    X <- qr.X(Q, ncol = length(Q$pivot))
    Z <- X[, seq_len(ncond), drop = FALSE]
    X <- X[, seq.int(ncond + 1L, ncol(X)), drop = FALSE]
    w <- attr(response, "RW")
    if (!is.null(w)) {
        X <- X / sqrt(w)
        Z <- Z / sqrt(w)
    }
    if (typed) {
        design <- attr(permutations, "control")
        strata <- getStrata(design, which = "plots")
        variables <- getPlots(design)$variables
        first <- which(!duplicated(strata))
        membership <- match(strata, strata[first])
        if (!isTRUE(all.equal(unname(X), unname(X[first[membership], , drop = FALSE]),
                              tolerance = 1e-10, check.attributes = FALSE)))
            stop("all constraints must be constant within each whole plot")
        profiles <- variables[first, , drop = FALSE]
        codes <- lapply(profiles, function(z) factor(match(z, unique(z))))
        groups <- do.call(interaction, c(codes, list(drop = TRUE)))
        for (group in levels(groups)) {
            take <- first[groups == group]
            expected <- X[rep(take[1L], length(take)), , drop = FALSE]
            if (!isTRUE(all.equal(unname(X[take, , drop = FALSE]), unname(expected),
                                  tolerance = 1e-10, check.attributes = FALSE)))
                stop("Plots(variables) must distinguish all distinct constraint rows")
        }
    }
    if (!is.null(x$pCCA)) {
        message("With conditioning variables, verify that the constraint permutation ",
                "scheme is justified conditional on those variables. ",
                "Condition() alone does not establish exchangeability.")
    }
    q <- x$CCA$qrank
    r <- N - Q$rank - 1L
    if (r < 1L)
        stop("the fitted model has no residual degrees of freedom")
    num <- x$CCA$tot.chi
    den <- x$CA$tot.chi
    F.0 <- (num / q) / (den / r)
    evaluate <- function(i) {
        constraintStatistic(response, X[indices[i, ], , drop = FALSE], Z)
    }
    nperm <- nrow(indices)
    if (is.null(parallel)) parallel <- 1L
    cluster <- inherits(parallel, "cluster")
    if (cluster) {
        values <- parLapply(parallel, seq_len(nperm), evaluate)
    } else if (parallel > 1L && .Platform$OS.type == "unix") {
        values <- mclapply(seq_len(nperm), evaluate, mc.cores = parallel)
    } else if (parallel > 1L) {
        workers <- makeCluster(parallel)
        on.exit(stopCluster(workers))
        values <- parLapply(workers, seq_len(nperm), evaluate)
    } else {
        values <- lapply(seq_len(nperm), evaluate)
    }
    values <- if (length(values)) do.call(rbind, values) else
        matrix(numeric(), nrow = 0L, ncol = 5L,
               dimnames = list(NULL, c("F", "num", "den", "q", "r")))
    sol <- list(call = match.call(), testcall = x$call, model = "constraints",
                F.0 = c(Model = F.0),
                F.perm = matrix(values[, 'F'], ncol = 1L),
                chi = c(Model = num, Residual = den),
                num = values[, 'num'], den = values[, 'den'],
                df = c(q, r), nperm = nperm,
                method = x$method, inertia = x$inertia, first = FALSE,
                termlabels = "Model", Random.seed = attr(permutations, "seed"),
                control = attr(permutations, "control"))
    ## An enumeration containing the observed arrangement already counts it.
    sol$observed <- isTRUE(attr(permutations, "observed"))
    control <- attr(permutations, "control")
    if (!is.null(control) && getObserved(control) &&
        !is.null(getAllperms(control)) && nperm == numPerms(N, control))
        sol$observed <- TRUE
    class(sol) <- "permutest.cca"
    sol
}
