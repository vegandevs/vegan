ordicomparisons <- function(object, ...)
    UseMethod("ordicomparisons")

ordicomparisons.default <- function(object, ...)
    stop("no 'ordicomparisons' method for objects of class ",
         paste(class(object), collapse = "/"))

ordicontrasts <- function(object, ...)
    UseMethod("ordicontrasts")

ordicontrasts.default <- function(object, ...)
    stop("no 'ordicontrasts' method for objects of class ",
         paste(class(object), collapse = "/"))

ordicontrasts.cca <- function(object, variable, newdata = NULL, by = NULL,
                                wts = NULL, ...)
{
    if (is.null(object$terms) || is.null(object$terminfo))
        stop("comparison profiles need a model fitted with the formula interface")
    if (length(variable) != 1L || !is.character(variable))
        stop("'variable' must name exactly one factor")
    dat <- .comparisonData(object)
    if (!variable %in% names(dat))
        stop(gettextf("variable '%s' was not found in the model data",
                      variable))
    if (!is.factor(dat[[variable]]))
        stop(gettextf("variable '%s' must be a factor", variable))
    predictors <- all.vars(.comparisonFormula(object))
    if (!variable %in% predictors)
        stop(gettextf("variable '%s' is not a constraint in the model",
                      variable))
    ## Formula symbols can also name settings such as a polynomial degree.
    predictors <- intersect(predictors, names(dat))
    if (!is.null(by)) {
        if (!is.character(by) || any(!by %in% names(dat)) || anyDuplicated(by))
            stop("'by' must name variables in the model data")
        if (variable %in% by)
            stop("the focal variable cannot also be a 'by' variable")
    }

    Q <- object$CCA$QR
    p <- length(object$terminfo$assign)
    if (is.null(Q) || !p)
        stop("model has no constraints")
    if (object$CCA$qrank < p)
        stop("comparison profiles are not available for aliased constraint matrices")
    ## Undo pivoting to recover the original constraint column names.
    nz <- length(Q$pivot) - p
    columns <- colnames(Q$qr)[order(Q$pivot)][nz + seq_len(p)]

    base <- .comparisonNewdata(predictors, dat, variable, by, newdata)
    split <- .comparisonSplit(base, by)
    meta <- .profileMetadata(predictors, dat, base, variable, by, newdata, wts)
    covariates <- intersect(intersect(names(base), names(dat)), predictors)
    covariates <- Filter(function(x) is.numeric(x) && any(!is.na(x)),
                         dat[covariates])
    ranges <- lapply(covariates, range, na.rm = TRUE)
    lev <- levels(dat[[variable]])
    ordered <- is.ordered(dat[[variable]])
    refs <- lapply(seq_along(split), function(g) {
        take <- split[[g]]
        gd <- base[take, , drop = FALSE]
        gw <- .comparisonWeights(gd, wts, take, nrow(base))
        R <- .comparisonRows(object, gd, variable, lev, ordered, gw, columns)
        group <- if (length(by)) gd[1L, by, drop = FALSE] else NULL
        extra <- .comparisonExtrapolation(gd, ranges)
        list(reference.rows = R, group = group, weights = gw,
             n.reference = nrow(gd), extrapolated = any(extra),
             extrapolation = extra)
    })

    labels <- make.unique(c("Profile", "N.reference", "Extrapolated",
                            variable, by))
    tabs <- lapply(refs, function(z) {
        tab <- data.frame(Profile = rep.int("", length(lev)),
                          stringsAsFactors = FALSE)
        tab[[labels[4L]]] <- factor(lev, levels = lev, ordered = ordered)
        if (!is.null(z$group)) {
            group <- z$group[rep.int(1L, length(lev)), , drop = FALSE]
            names(group) <- labels[-seq_len(4L)]
            tab <- cbind(tab, group, row.names = NULL)
        }
        tab$N.reference <- rep.int(z$n.reference, length(lev))
        tab$Extrapolated <- rep.int(z$extrapolated, length(lev))
        tab
    })
    out <- do.call(rbind, tabs)
    rownames(out) <- NULL
    ids <- paste0("b", seq_len(nrow(out)))
    clash <- ids %in% lev
    if (any(clash)) {
        ids[clash] <- paste0(".", ids[clash])
        ids <- make.unique(c(lev, ids))[length(lev) + seq_along(ids)]
    }
    out$Profile <- ids
    R <- do.call(rbind, lapply(refs, `[[`, "reference.rows"))
    rownames(R) <- out$Profile
    group.index <- rep(seq_along(refs), each = length(lev))

    ## Profile identifiers link the displayed rows to the fitted design rows.
    attr(out, "profile.matrix") <- R
    attr(out, "group.index") <- group.index
    attr(out, "reference") <- refs
    attr(out, "reference.data") <- base
    attr(out, "model") <- object
    attr(out, "variable") <- variable
    attr(out, "levels") <- lev
    attr(out, "by") <- by
    attr(out, "newdata") <- newdata
    attr(out, "weights") <- wts
    attr(out, "population") <- meta$population
    attr(out, "averaged") <- meta$averaged
    attr(out, "fixed") <- meta$fixed
    attr(out, "weight.label") <- meta$weights
    attr(out, "extrapolation") <- .comparisonExtrapolation(base, ranges)
    class(out) <- c("vegan_contrasts", "data.frame")
    out
}

print.vegan_contrasts <- function(x, ...)
{
    none <- function(z) if (length(z)) paste(z, collapse = ", ") else "none"
    cat("Comparison profiles for constrained ordination\n\n")
    cat("Focal variable:      ", attr(x, "variable"), "\n")
    if (length(attr(x, "by")))
        cat("By:                  ", paste(attr(x, "by"), collapse = ", "),
            "\n")
    cat("Reference population:", attr(x, "population"), "\n")
    cat("Averaged over:       ", none(attr(x, "averaged")), "\n")
    if (length(attr(x, "fixed")))
        cat("Fixed at:            ", none(attr(x, "fixed")), "\n")
    cat("Weights:             ", attr(x, "weight.label"), "\n\n")
    print.data.frame(x, row.names = FALSE, ...)
    invisible(x)
}

ordicomparisons.cca <-
    function(object, variable, comparison = "pairwise", reference = NULL,
             newdata = NULL, by = NULL, joint = FALSE,
             permutations = how(nperm = 999), strata = NULL,
             parallel = getOption("mc.cores"), adjust = "none", wts = NULL,
             ...)
{
    profiles <- ordicontrasts(object, variable = variable,
                                 newdata = newdata, by = by, wts = wts)
    ordicomparisons(profiles, comparison = comparison, reference = reference,
                joint = joint, permutations = permutations, strata = strata,
                parallel = parallel, adjust = adjust, ...)
}

ordicomparisons.vegan_contrasts <-
    function(object, comparison = "pairwise", reference = NULL,
             joint = FALSE, permutations = how(nperm = 999), strata = NULL,
             parallel = getOption("mc.cores"), adjust = "none", ...)
{
    adjust <- match.arg(adjust, p.adjust.methods)
    if (length(joint) != 1L || is.na(joint))
        stop("'joint' must be TRUE or FALSE")

    fit <- attr(object, "model")
    R <- attr(object, "profile.matrix")
    ids <- object$Profile
    lev <- attr(object, "levels")
    by <- attr(object, "by")
    groups <- attr(object, "group.index")
    if (is.null(R) || nrow(R) != nrow(object) ||
        !identical(rownames(R), ids) || length(groups) != nrow(object))
        stop("profile rows no longer match their retained model matrix; ",
             "reconstruct them with 'ordicontrasts()'")
    Q <- fit$CCA$QR
    p <- ncol(R)
    nz <- length(Q$pivot) - p
    tested <- Q$rank - p + seq_len(p)
    pivot <- Q$pivot[tested] - nz
    triangular <- qr.R(Q)[tested, tested, drop = FALSE]
    basis <- qr.qy(Q, diag(nrow = nrow(Q$qr), ncol = Q$rank))
    n <- nrow(basis)
    if (is.matrix(permutations) &&
        (!is.numeric(permutations) || is.complex(permutations) ||
         ncol(permutations) != n ||
         anyNA(permutations) || any(permutations < 1L | permutations > n) ||
         any(apply(permutations, 1L, anyDuplicated) > 0L)))
        stop("each row of 'permutations' must contain each observation once")
    perm <- getPermuteMatrix(permutations, n, strata = strata)

    global <- .comparisonUsesProfileRows(comparison, ids, length(lev))
    if (global) {
        A <- .comparisonContrasts(comparison, ids, reference = NULL)
        L <- .comparisonMatrix(A, R)
        specs <- lapply(seq_len(nrow(A)), function(i)
            list(label = rownames(A)[i], L = L[i, , drop = FALSE],
                 group = NULL))
    } else {
        specs <- unlist(lapply(unique(groups), function(g) {
            take <- which(groups == g)
            A <- .comparisonContrasts(comparison, lev, reference)
            L <- .comparisonMatrix(A, R[take, , drop = FALSE])
            group <- attr(object, "reference")[[g]]$group
            lapply(seq_len(nrow(A)), function(i)
                list(label = rownames(A)[i],
                     L = L[i, , drop = FALSE], group = group))
        }), recursive = FALSE)
    }
    if (joint && length(specs) > 1L) {
        L <- do.call(rbind, lapply(specs, `[[`, "L"))
        specs <- list(list(label = paste(vapply(specs, `[[`, "", "label"),
                                             collapse = "; "),
                           L = L, group = NULL))
    }

    ans <- lapply(specs, function(spec) {
        z <- .comparisonTest(fit, basis, triangular, pivot, spec$L, perm,
                             parallel)
        tab <- data.frame(Comparison = spec$label, Df = z$df,
                          Inertia = z$inertia,
                          R2 = z$inertia / fit$tot.chi,
                          F = z$F, N.Perm = nrow(perm),
                          check.names = FALSE)
        tab[["Pr(>F)"]] <- z$p
        if (!is.null(group <- spec$group)) {
            names(group) <- make.unique(c(names(tab), names(group)))[-seq_along(tab)]
            tab <- cbind(group, tab, row.names = NULL)
        }
        tab
    })
    out <- do.call(rbind, ans)
    rownames(out) <- NULL
    if (adjust != "none")
        out[["Pr(>F).adj"]] <- p.adjust(out[["Pr(>F)"]], method = adjust)
    attr(out, "comparison.matrices") <- list(
        profiles = R, comparisons = lapply(specs, `[[`, "L"))
    attr(out, "contrasts") <- object
    attr(out, "model") <- fit
    attr(out, "estimable") <- rep.int(TRUE, nrow(out))
    attr(out, "extrapolation") <- attr(object, "extrapolation")
    attr(out, "variable") <- attr(object, "variable")
    attr(out, "comparison") <- comparison
    attr(out, "newdata") <- attr(object, "newdata")
    attr(out, "by") <- by
    attr(out, "adjust") <- adjust
    attr(out, "weights") <- attr(object, "weights")
    attr(out, "joint") <- joint
    attr(out, "Random.seed") <- attr(perm, "seed")
    attr(out, "control") <- attr(perm, "control")
    class(out) <- c("vegan_comparisons", "data.frame")
    out
}

print.vegan_comparisons <- function(x, ...)
{
    cat("Permutation comparisons for constrained ordination\n")
    cat("Variable:", attr(x, "variable"), "\n\n")
    print.data.frame(x, ...)
    invisible(x)
}

.comparisonData <- function(object)
{
    call <- object$call
    if (is.null(call$data))
        stop("no sufficient information to reconstruct model data")
    dat <- eval(call$data, environment(object$terms), .GlobalEnv)
    if (inherits(dat, "environment"))
        dat <- as.data.frame(as.list(dat))
    else
        dat <- as.data.frame(dat)
    if (!is.null(object$subset))
        dat <- dat[object$subset, , drop = FALSE]
    if (!is.null(object$na.action))
        dat <- dat[-as.integer(object$na.action), , drop = FALSE]
    if (nrow(dat) != nrow(ordiYbar(object, "initial")))
        stop("could not align reconstructed model data with fitted ",
             "observations")
    xlev <- object$terminfo$xlev
    for (nm in intersect(names(xlev), names(dat)))
        dat[[nm]] <- factor(dat[[nm]], levels = xlev[[nm]],
                            ordered = isTRUE(object$terminfo$ordered[[nm]]))
    dat
}

.comparisonNewdata <- function(predictors, dat, variable, by, newdata)
{
    if (is.null(newdata))
        return(dat)
    if (is.character(newdata)) {
        choices <- c("balanced", "mean")
        if (length(newdata) != 1L || !newdata %in% choices)
            stop("character 'newdata' must be \"balanced\" or \"mean\"")
        vars <- unique(c(setdiff(predictors, variable), by))
        if (newdata == "mean")
            return(.comparisonMeanData(dat, vars, by))
        vals <- lapply(vars, function(v) {
            x <- dat[[v]]
            if (is.factor(x)) levels(x)
            else if (is.numeric(x)) mean(x)
            else unique(x)
        })
        names(vals) <- vars
        if (!length(vals))
            return(dat[1L, vars, drop = FALSE])
        return(expand.grid(vals, KEEP.OUT.ATTRS = FALSE,
                           stringsAsFactors = FALSE))
    }
    if (!is.data.frame(newdata))
        stop("'newdata' must be NULL, \"balanced\", \"mean\", or a data frame")
    newdata
}

.comparisonTypicalValue <- function(x)
{
    ok <- !is.na(x)
    if (!any(ok))
        stop("cannot construct a mean-or-mode profile from all-missing values")
    if (is.numeric(x))
        return(mean(x[ok]))
    if (is.factor(x)) {
        count <- tabulate(as.integer(x[ok]), nbins = nlevels(x))
        return(factor(levels(x)[which.max(count)], levels = levels(x),
                      ordered = is.ordered(x)))
    }
    observed <- x[ok]
    values <- unique(observed)
    values[which.max(tabulate(match(observed, values)))]
}

.comparisonMeanData <- function(dat, vars, by)
{
    if (!length(vars))
        return(dat[1L, vars, drop = FALSE])
    groups <- .comparisonSplit(dat, by)
    rows <- lapply(groups, function(take) {
        gd <- dat[take, , drop = FALSE]
        values <- lapply(vars, function(v) {
            if (v %in% by) gd[[v]][1L]
            else .comparisonTypicalValue(gd[[v]])
        })
        names(values) <- vars
        as.data.frame(values, optional = TRUE, stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
}

.profileMetadata <- function(predictors, dat, base, variable, by, newdata, wts)
{
    vars <- setdiff(predictors, c(variable, by))
    fixed <- character()
    if (is.null(newdata)) {
        population <- paste0("fitted analysis data (", nrow(base), " rows)")
        averaged <- vars
    } else if (is.character(newdata)) {
        if (newdata == "balanced") {
            population <- paste0("balanced reference grid (", nrow(base),
                                 " rows)")
            numeric <- vars[vapply(dat[vars], is.numeric, logical(1))]
            averaged <- setdiff(vars, numeric)
            if (length(numeric))
                fixed <- paste0(numeric, " = ",
                    vapply(base[numeric], function(x)
                        format(signif(x[[1L]], 6L), trim = TRUE),
                        character(1)))
        } else {
            what <- if (nrow(base) == 1L) "profile" else "profiles"
            population <- paste0("mean-or-mode reference ", what, " (",
                                 nrow(base), if (nrow(base) == 1L) " row)"
                                 else " rows)")
            averaged <- character()
            if (length(vars)) {
                if (length(by))
                    fixed <- paste0(vars, " = mean/mode within ",
                                    paste(by, collapse = ", "))
                else
                    fixed <- paste0(vars, " = ", vapply(base[vars],
                        function(x) {
                            z <- x[[1L]]
                            if (is.numeric(z))
                                format(signif(z, 6L), trim = TRUE)
                            else as.character(z)
                        }, character(1)))
            }
        }
    } else {
        population <- paste0("user-supplied data (", nrow(base), " rows)")
        present <- intersect(vars, names(base))
        is.fixed <- vapply(base[present], function(x)
            length(unique(x[!is.na(x)])) <= 1L, logical(1))
        fixed.vars <- present[is.fixed]
        averaged <- setdiff(vars, fixed.vars)
        if (length(fixed.vars))
            fixed <- paste0(fixed.vars, " = ", vapply(base[fixed.vars],
                function(x) as.character(x[which(!is.na(x))[1L]]),
                character(1)))
    }
    if (length(by) && length(averaged))
        averaged <- paste0(averaged, " within ", paste(by, collapse = ", "))
    weights <- if (is.null(wts)) "equal" else if (is.character(wts)) wts
               else "supplied numeric weights"
    if (length(by))
        weights <- paste(weights, "within", paste(by, collapse = ", "),
                         "groups")
    list(population = population, averaged = averaged, fixed = fixed,
         weights = weights)
}

.comparisonExtrapolation <- function(newdata, ranges)
{
    vapply(names(ranges), function(v) {
        r <- ranges[[v]]
        any(newdata[[v]] < r[1L] | newdata[[v]] > r[2L], na.rm = TRUE)
    }, logical(1))
}

.comparisonSplit <- function(dat, by)
{
    if (!nrow(dat))
        stop("'newdata' must contain at least one row")
    if (!length(by))
        return(list(seq_len(nrow(dat))))
    if (any(!by %in% names(dat)))
        stop("the reference population must contain all 'by' variables")
    if (anyNA(dat[by], recursive = TRUE))
        stop("'by' variables must not contain missing values")
    ## Integer codes keep different groups with the same pasted label apart.
    group <- lapply(dat[by], function(x) as.integer(as.factor(x)))
    split(seq_len(nrow(dat)), interaction(group, drop = TRUE, lex.order = TRUE))
}

.comparisonWeights <- function(dat, wts, take, n)
{
    if (is.null(wts))
        w <- rep.int(1, nrow(dat))
    else if (is.character(wts) && length(wts) == 1L && wts %in% names(dat))
        w <- dat[[wts]]
    else if (is.numeric(wts) && length(wts) == n)
        w <- wts[take]
    else if (is.numeric(wts) && length(wts) == nrow(dat))
        w <- wts
    else
        stop("'wts' must be a numeric vector or a column in 'newdata'")
    if (!is.numeric(w) || any(!is.finite(w)) || any(w < 0) || !any(w > 0))
        stop("comparison weights must be non-negative, finite, and sum ",
             "above zero")
    w <- w / max(w)
    w / sum(w)
}

.comparisonFormula <- function(object)
{
    Terms <- terms(object)
    labs <- attr(Terms, "term.labels")
    labs <- labs[!grepl("^Condition\\(", labs)]
    reformulate(labs, intercept = attr(Terms, "intercept"),
                env = environment(Terms))
}

.comparisonRows <- function(object, dat, variable, levels, ordered, w, columns)
{
    Terms <- object$terminfo$constraint.terms
    if (is.null(Terms)) {
        f <- .comparisonFormula(object)
        vars <- attr(terms(f), "variables")[-1L]
        if (any(vapply(vars, is.call, logical(1))))
            stop("refit the model to retain fitted transformations for comparisons")
        Terms <- terms(f)
    }
    ctr <- object$terminfo$contrasts
    xlev <- object$terminfo$xlev
    variables <- names(attr(Terms, "dataClasses"))
    if (is.null(variables))
        variables <- all.vars(Terms)
    xlev <- xlev[names(xlev) %in% variables]
    if (length(xlev) && is.null(ctr))
        stop("refit the model to retain its factor contrast coding for comparisons")
    rows <- lapply(levels, function(lev) {
        nd <- dat
        nd[[variable]] <- factor(lev, levels = levels,
                                 ordered = ordered)
        mf <- model.frame(Terms, nd, na.action = na.fail, xlev = xlev)
        mm <- model.matrix(Terms, mf, contrasts.arg = ctr)
        mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
        miss <- setdiff(columns, colnames(mm))
        if (length(miss))
            stop("could not reconstruct fitted model-matrix columns: ",
                 paste(miss, collapse = ", "))
        drop(crossprod(w, mm[, columns, drop = FALSE]))
    })
    out <- do.call(rbind, rows)
    rownames(out) <- levels
    colnames(out) <- columns
    out
}

## Explicit profile identifiers can cross groups. Factor-level recipes and
## equations are evaluated separately within each group.
.comparisonUsesProfileRows <- function(comparison, ids, nlevels)
{
    recipes <- c("pairwise", "sequential", "consecutive", "reference",
                 "trt.vs.ctrl", "deviation", "effect", "helmert", "poly")
    if (is.character(comparison)) {
        if (length(comparison) == 1L && comparison %in% recipes)
            return(FALSE)
        return(any(vapply(comparison, function(z) {
            expr <- try(str2lang(z), silent = TRUE)
            !inherits(expr, "try-error") && any(all.vars(expr) %in% ids)
        }, logical(1))))
    }
    if (is.list(comparison))
        return(any(vapply(comparison, .comparisonUsesProfileRows, logical(1),
                          ids = ids, nlevels = nlevels)))
    if (is.numeric(comparison)) {
        n <- if (is.matrix(comparison)) ncol(comparison) else length(comparison)
        keys <- if (is.matrix(comparison)) colnames(comparison)
                else names(comparison)
        if (!is.null(keys) && any(keys %in% ids))
            return(TRUE)
        return(length(ids) != nlevels && n == length(ids))
    }
    FALSE
}

.comparisonContrasts <- function(comparison, keys, reference = NULL)
{
    k <- length(keys)
    recipe <- length(comparison) == 1L && is.character(comparison) &&
        comparison %in% c("pairwise", "sequential", "consecutive",
                          "reference", "trt.vs.ctrl", "deviation",
                          "effect", "helmert", "poly")
    if (recipe) {
        if (k < 2L)
            stop("at least two quantities are needed for comparisons")
        if (comparison == "pairwise") {
            cmb <- combn(k, 2L)
            A <- matrix(0, ncol = k, nrow = ncol(cmb))
            for (i in seq_len(nrow(A))) A[i, cmb[, i]] <- c(-1, 1)
        } else if (comparison %in% c("sequential", "consecutive")) {
            A <- matrix(0, ncol = k, nrow = k - 1L)
            for (i in seq_len(k - 1L)) A[i, i:(i + 1L)] <- c(-1, 1)
        } else if (comparison %in% c("reference", "trt.vs.ctrl")) {
            if (is.null(reference)) reference <- keys[[1L]]
            ref <- match(reference, keys)
            if (length(ref) != 1L || is.na(ref))
                stop("'reference' must name one of the quantities")
            take <- setdiff(seq_len(k), ref)
            A <- matrix(0, nrow = length(take), ncol = k)
            for (i in seq_along(take)) A[i, c(ref, take[i])] <- c(-1, 1)
        } else if (comparison %in% c("deviation", "effect")) {
            A <- diag(k) - matrix(1 / k, k, k)
        } else if (comparison == "helmert") {
            A <- t(stats::contr.helmert(k))
        } else {
            A <- t(stats::contr.poly(k))
        }
        colnames(A) <- keys
        rownames(A) <- .comparisonLabels(A, keys)
        return(A)
    }
    if (is.character(comparison)) {
        A <- do.call(rbind, lapply(comparison, .comparisonExpression,
                                   keys = keys))
        rownames(A) <- if (!is.null(names(comparison)) &&
                           all(nzchar(names(comparison)))) names(comparison)
                      else comparison
    } else if (is.list(comparison)) {
        matrices <- lapply(comparison, .comparisonContrasts, keys = keys,
                           reference = reference)
        labels <- names(matrices)
        for (i in which(!is.na(labels) & nzchar(labels))) {
            mat <- matrices[[i]]
            rownames(mat) <- if (nrow(mat) == 1L) labels[i]
                             else paste(labels[i], rownames(mat), sep = ": ")
            matrices[[i]] <- mat
        }
        A <- do.call(rbind, unname(matrices))
    } else if (is.numeric(comparison)) {
        if (is.matrix(comparison)) {
            A <- comparison
        } else {
            A <- matrix(comparison, nrow = 1L,
                        dimnames = list(NULL, names(comparison)))
        }
        if (!is.null(colnames(A))) {
            if (ncol(A) != k || anyDuplicated(colnames(A)) ||
                !setequal(colnames(A), keys))
                stop("named contrast columns must match the quantities exactly")
            A <- A[, keys, drop = FALSE]
        } else if (ncol(A) != k) {
            stop("numeric comparisons must have one column per quantity")
        }
        if (is.null(rownames(A))) rownames(A) <- .comparisonLabels(A, keys)
    } else {
        stop("unsupported 'comparison' specification")
    }
    colnames(A) <- keys
    A
}

.comparisonMatrix <- function(A, R)
{
    if (!length(A) || any(!is.finite(A)))
        stop("comparison weights must be finite and non-empty")
    L <- matrix(0, nrow(A), ncol(R), dimnames = list(rownames(A), colnames(R)))
    for (i in seq_len(nrow(A))) {
        a <- A[i, ]
        take <- which(a != 0)
        if (!length(take))
            next
        a <- a[take]
        anchor <- which.max(abs(a))
        a <- a / abs(a[anchor])
        if (abs(sum(a)) > 8 * .Machine$double.eps * sum(abs(a)))
            stop("comparison weights must sum to zero")
        ## Subtract a common profile before multiplication so a large
        ## shared offset cannot swamp a small difference.
        diffs <- sweep(R[take, , drop = FALSE], 2L, R[take[anchor], ], "-")
        if (any(!is.finite(diffs)))
            stop("comparison produced non-finite model-matrix differences")
        x <- drop(crossprod(a, diffs))
        if (any(!is.finite(x)))
            stop("comparison produced non-finite model-matrix differences")
        ## Ignore cancellation at the rounding precision of the product.
        bound <- drop(crossprod(8 * .Machine$double.eps * abs(a), abs(diffs)))
        x[abs(x) <= bound] <- 0
        L[i, ] <- x
    }
    L
}

.comparisonExpression <- function(expression, keys)
{
    if (length(expression) != 1L)
        stop("each comparison equation must be one character string")
    expr <- tryCatch(str2lang(expression), error = function(e)
        stop("invalid comparison expression: ", conditionMessage(e)))
    zero <- setNames(numeric(length(keys)), keys)
    walk <- function(e) {
        if (is.numeric(e) && length(e) == 1L)
            return(list(coef = zero, constant = as.numeric(e)))
        if (is.name(e)) {
            nm <- as.character(e)
            if (!nm %in% keys)
                stop(gettextf("unknown comparison name '%s'", nm))
            z <- zero
            z[[nm]] <- 1
            return(list(coef = z, constant = 0))
        }
        if (!is.call(e)) stop("comparison expressions must be linear")
        op <- as.character(e[[1L]])
        if (op == "(")
            return(walk(e[[2L]]))
        if (!op %in% c("+", "-", "*", "/", "="))
            stop(gettextf("operator '%s' is not allowed", op))
        if (length(e) == 2L && op %in% c("+", "-")) {
            a <- walk(e[[2L]])
            if (op == "-") {
                a$coef <- -a$coef
                a$constant <- -a$constant
            }
            return(a)
        }
        if (length(e) != 3L) stop("invalid linear comparison expression")
        a <- walk(e[[2L]])
        b <- walk(e[[3L]])
        if (any(!is.finite(c(a$coef, a$constant, b$coef, b$constant))))
            stop("comparison expressions must have finite weights and constants")
        if (op %in% c("+", "-", "=")) {
            s <- if (op == "+") 1 else -1
            return(list(coef = a$coef + s * b$coef,
                        constant = a$constant + s * b$constant))
        }
        if (op == "*") {
            if (all(a$coef == 0))
                return(list(coef = a$constant * b$coef,
                            constant = a$constant * b$constant))
            if (all(b$coef == 0))
                return(list(coef = b$constant * a$coef,
                            constant = b$constant * a$constant))
            stop("multiplication is only allowed by numeric constants")
        }
        if (any(b$coef != 0) || b$constant == 0)
            stop("division is only allowed by a non-zero numeric constant")
        list(coef = a$coef / b$constant,
             constant = a$constant / b$constant)
    }
    z <- walk(expr)
    if (any(!is.finite(c(z$coef, z$constant))))
        stop("comparison expressions must have finite weights and constants")
    if (z$constant != 0)
        stop("only comparisons with a zero right-hand side are supported")
    unname(z$coef)
}

.comparisonLabels <- function(A, keys)
{
    apply(A, 1L, function(a) {
        pos <- which(a > 0)
        neg <- which(a < 0)
        fmt <- function(i) paste0(ifelse(abs(a[i]) == 1, "",
                                         paste0(abs(a[i]), "*")),
                                  keys[i], collapse = " + ")
        paste0(fmt(pos), " - ", fmt(neg))
    })
}

.comparisonTest <- function(object, basis, triangular, pivot, L, permutations,
                            parallel)
{
    if (!any(L != 0))
        stop("comparison has rank zero")
    if (nrow(basis) - ncol(basis) - 1L < 1L)
        stop("comparison test has no residual degrees of freedom")

    ## If the residualised constraints are Qx R, the hypothesis in Qx
    ## coordinates is L R^-1. Use a triangular solve, respecting pivoting.
    U <- backsolve(triangular, t(L[, pivot, drop = FALSE]), transpose = TRUE)
    if (any(!is.finite(U)))
        stop("comparison is not estimable at the scale of the fitted design")
    size <- apply(abs(U), 2L, max)
    nonzero <- size > 0
    if (!any(nonzero))
        stop("comparison has rank zero")
    U <- sweep(U[, nonzero, drop = FALSE], 2L, size[nonzero], "/")
    rotation <- qr(U, tol = sqrt(.Machine$double.eps))
    rank <- rotation$rank
    constraints <- ncol(basis) - nrow(triangular) + seq_len(nrow(triangular))
    rotated <- basis[, constraints, drop = FALSE] %*%
        qr.Q(rotation, complete = TRUE)
    tested <- seq_len(rank)
    Z0 <- cbind(basis[, -constraints, drop = FALSE],
                rotated[, -tested, drop = FALSE])
    Q <- qr(cbind(Z0, rotated[, tested, drop = FALSE]))
    QZ <- if (ncol(Z0)) qr(Z0) else NULL

    Y <- ordiYbar(object, "initial")
    qty <- qr.qty(Q, Y)
    tol <- nrow(Y) * .Machine$double.eps
    if (inherits(object, "dbrda")) {
        inertia <- diag(qr.qty(Q, t(qty)))
        min.residual <- tol * norm(Y, "F")
    } else {
        inertia <- rowSums(qty^2)
        min.residual <- tol^2 * sum(inertia)
    }
    contrast <- sum(inertia[ncol(Z0) + tested])
    residual <- sum(inertia[seq.int(Q$rank + 1L, nrow(Y))])
    ## A perfectly fitted response can leave a round-off-sized residual.
    if (!is.finite(residual) || abs(residual) <= min.residual)
        stop("comparison test has zero or non-finite residual inertia")

    ## Only the design QR and inertia are needed by permutest.cca().
    ## The ranks describe design spaces regardless of the response inertia.
    fit <- list(Ybar = Y, call = object$call, method = object$method,
                inertia = object$inertia, comparison = TRUE,
                CCA = list(rank = rank, qrank = rank, QR = Q, tot.chi = contrast),
                CA = list(tot.chi = residual))
    if (!is.null(QZ))
        fit$pCCA <- list(rank = QZ$rank, QR = QZ)
    class(fit) <- class(object)
    tst <- permutest.cca(fit, permutations = permutations,
                         model = "reduced", parallel = parallel)
    eps <- sqrt(.Machine$double.eps)
    p <- (sum(tst$F.perm >= tst$F.0 - eps) + 1) / (tst$nperm + 1)
    list(df = rank, inertia = contrast, F = unname(tst$F.0), p = p)
}

.comparisonsAnova <- function(x, object)
{
    columns <- c("Df", "Inertia", "F", "Pr(>F)")
    if ("Pr(>F).adj" %in% names(x))
        columns <- c(columns, "Pr(>F).adj")
    keep <- x[, columns, drop = FALSE]
    varname <- if (inherits(object, c("capscale", "dbrda")) &&
                   object$adjust == 1)
        "SumOfSqs" else if (inherits(object, "rda")) "Variance" else "ChiSquare"
    names(keep)[2L] <- varname
    labels <- x$Comparison
    group <- seq_len(match("Comparison", names(x)) - 1L)
    if (length(group)) {
        values <- Map(function(nm, value) paste(nm, value, sep = " = "),
                      attr(x, "by"), x[group])
        labels <- paste0(labels, " (",
                         do.call(paste, c(unname(values), sep = ", ")), ")")
    }
    rownames(keep) <- make.unique(labels)
    attr(keep, "heading") <- c("Permutation test for specified comparisons\n",
                               paste("Model:", c(object$call)))
    attr(keep, "Random.seed") <- attr(x, "Random.seed")
    attr(keep, "control") <- attr(x, "control")
    class(keep) <- c("anova.cca", "anova", "data.frame")
    keep
}
