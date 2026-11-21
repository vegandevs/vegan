## Align the design to the fitted predictors before generating assignments.
constraintPermutationDesign <- function(x, permutations)
{
    if (!inherits(permutations, "how"))
        return(permutations)
    plots <- getPlots(permutations)
    data <- plots$variables
    if (is.null(data))
        return(permutations)
    qr <- x$CCA$QR$qr
    target <- rownames(qr)
    if (!is.null(target) && .row_names_info(data, 1L) > 0L) {
        ## The response can have different row names from its predictor data.
        take <- match(target, rownames(data))
    } else {
        take <- seq_len(nrow(data))
        if (length(take) != nrow(qr)) {
            if (!is.null(x$subset))
                take <- take[x$subset]
            if (!is.null(x$na.action))
                take <- take[-x$na.action]
        }
    }
    if (length(take) != nrow(qr) || anyNA(take))
        stop("plot variables cannot be aligned with the fitted observations")
    if (identical(take, seq_len(nrow(data))))
        return(permutations)
    plots <- update(plots, strata = getStrata(plots)[take],
                    variables = data[take, , drop = FALSE])
    blocks <- getBlocks(permutations)
    if (!is.null(blocks))
        blocks <- blocks[take]
    ## A cached enumeration was generated for the original set of rows.
    update(permutations, plots = plots, blocks = blocks, all.perms = NULL)
}

## Use the same eigenvalue cutoff as ordConstrained(). adonis0() uses
## the raw trace instead.
constraintInertia <- function(Y, isDB, isAdonis)
{
    if (isAdonis)
        return(sum(diag(Y)))
    eig <- if (isDB) eigen(Y, symmetric = TRUE, only.values = TRUE)$values else
        svd(Y, nu = 0L, nv = 0L)$d^2
    zero <- sqrt(.Machine$double.eps)
    sum(eig[abs(eig) >= max(zero, zero * eig[1L])])
}

constraintStatistic <- function(Y, X, Z, w, isDB, isAdonis)
{
    if (is.null(w)) {
        X <- scale(X, scale = FALSE)
    } else {
        centre <- apply(X, 2, weighted.mean, w = w)
        X <- sqrt(w) * scale(X, center = centre, scale = FALSE)
    }
    Q <- qr(cbind(Z, X))
    kept <- seq_len(Q$rank)
    r <- nrow(Y) - Q$rank - 1L
    if (r < 1L)
        stop("a constraint permutation leaves no residual degrees of freedom")
    fit <- kept[Q$pivot[kept] > ncol(Z)]
    res <- seq.int(Q$rank + 1L, nrow(Y))
    q <- length(fit)
    ## In QR coordinates, constraints and residuals occupy disjoint
    ## rows (and columns for Gower-centred dissimilarities).
    Y <- qr.qty(Q, Y)
    if (isDB) {
        Y <- qr.qty(Q, t(Y))
        fitted <- Y[fit, fit, drop = FALSE]
        residual <- Y[res, res, drop = FALSE]
    } else {
        fitted <- Y[fit, , drop = FALSE]
        residual <- Y[res, , drop = FALSE]
    }
    num <- if (q > 0L) constraintInertia(fitted, isDB, isAdonis) else 0
    den <- constraintInertia(residual, isDB, isAdonis)
    F <- if (q > 0L) (num / q) / (den / r) else 0
    c(F = F, num = num, den = den)
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
    if (!is.null(w))
        X <- X / sqrt(w)
    if (typed) {
        design <- attr(permutations, "control")
        plot <- getStrata(design, which = "plots")
        variables <- getPlots(design)$variables
        target <- rownames(Q$qr)
        if (!is.null(target) && .row_names_info(variables, 1L) > 0L &&
            !identical(rownames(variables), target))
            stop("plot permutations must use the fitted order of observations")
        first <- which(!duplicated(plot))
        membership <- match(plot, plot[first])
        ## Allow QR reconstruction error on the scale of each column.
        tol <- 1e-10 * apply(abs(X), 2L, max)
        delta <- X - X[first[membership], , drop = FALSE]
        if (any(sweep(abs(delta), 2L, tol, ">")))
            stop("all constraints must be constant within each whole plot")
        ## Plots with the same assignment profile must be interchangeable.
        codes <- lapply(variables[first, , drop = FALSE],
                        function(z) match(z, unique(z)))
        profile <- interaction(codes, drop = TRUE)
        delta <- X[first, , drop = FALSE] -
            X[first[match(profile, profile)], , drop = FALSE]
        if (any(sweep(abs(delta), 2L, tol, ">")))
            stop("Plots(variables) must distinguish all distinct constraint rows")
    }
    if (!is.null(x$pCCA)) {
        message("With conditioning variables, verify that the constraint permutation ",
                "scheme is justified conditional on those variables. ",
                "Condition() alone does not establish exchangeability.")
    }
    isDB <- inherits(x, "dbrda")
    isAdonis <- inherits(x, "adonis2")
    if (ncond > 0L && QZ$rank > 0L) {
        response <- qr.resid(QZ, response)
        if (isDB)
            response <- qr.resid(QZ, t(response))
    }
    q <- x$CCA$qrank
    r <- N - Q$rank - 1L
    if (r < 1L)
        stop("the fitted model has no residual degrees of freedom")
    num <- x$CCA$tot.chi
    den <- x$CA$tot.chi
    F.0 <- (num / q) / (den / r)
    evaluate <- function(i) {
        constraintStatistic(response, X[indices[i, ], , drop = FALSE],
                            Z, w, isDB, isAdonis)
    }
    nperm <- nrow(indices)
    if (is.null(parallel)) parallel <- 1L
    cluster <- inherits(parallel, "cluster")
    if (cluster) {
        values <- parLapply(parallel, seq_len(nperm), evaluate)
    } else if (parallel > 1L && .Platform$OS.type == "unix") {
        evaluateFork <- function(i) tryCatch(evaluate(i), error = identity)
        values <- mclapply(seq_len(nperm), evaluateFork, mc.cores = parallel)
        failed <- which(vapply(values, inherits, logical(1), "error"))
        if (length(failed))
            stop(values[[failed[1L]]])
        if (any(vapply(values, is.null, logical(1))))
            stop("a constraint permutation worker returned no result")
    } else if (parallel > 1L) {
        workers <- makeCluster(parallel)
        on.exit(stopCluster(workers))
        values <- parLapply(workers, seq_len(nperm), evaluate)
    } else {
        values <- lapply(seq_len(nperm), evaluate)
    }
    values <- if (nperm > 0L) do.call(rbind, values) else
        matrix(numeric(), nrow = 0L, ncol = 3L,
               dimnames = list(NULL, c("F", "num", "den")))
    control <- attr(permutations, "control")
    ## An enumeration containing the observed arrangement already counts it.
    observed <- isTRUE(attr(permutations, "observed"))
    if (!observed && inherits(control, "how") && getObserved(control) &&
        !is.null(getAllperms(control)) && nperm == numPerms(N, control))
        observed <- TRUE
    sol <- list(call = match.call(), testcall = x$call, model = "constraints",
                F.0 = c(Model = F.0), F.perm = matrix(values[, "F"], ncol = 1L),
                chi = c(Model = num, Residual = den),
                num = values[, "num"], den = values[, "den"],
                df = c(q, r), nperm = nperm, method = x$method,
                inertia = x$inertia, first = FALSE, termlabels = "Model",
                Random.seed = attr(permutations, "seed"), control = control,
                observed = observed)
    class(sol) <- "permutest.cca"
    sol
}
