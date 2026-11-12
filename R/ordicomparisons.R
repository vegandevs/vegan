#' Comparisons for constrained ordination models
#'
#' Construct scientifically meaningful comparisons of factor levels and test
#' them with contrast-specific permutation pseudo-F statistics.
#'
#' @param object A fitted constrained ordination model or a profile table
#'   returned by `ordicontrasts()`.
#' @param variable The name of one focal factor in the model.
#' @param comparison A built-in comparison name, character equations, or a
#'   numeric contrast vector or matrix.
#' @param reference Optional reference or control level. The first fitted
#'   factor level is used when this is omitted.
#' @param joint Test all resulting comparisons jointly.
#' @param permutations A permutation count, matrix, or `permute::how()` object.
#' @param strata Optional permutation strata.
#' @param parallel Number of parallel workers or a socket cluster.
#' @param adjust A method accepted by `stats::p.adjust()`.
#' @param ... Arguments passed to methods.
#'
#' @return A data frame of class `vegan_comparisons`. The `contrasts`
#'   attribute retains the `vegan_contrasts` object used to define the tests.
#' @export
ordicomparisons <- function(object, ...)
    UseMethod("ordicomparisons")

#' Default comparisons method
#'
#' Report that no comparisons method is available for an object.
#'
#' @param object An unsupported object.
#' @param ... Unused arguments.
#'
#' @return This method always throws an error.
#' @keywords internal
ordicomparisons.default <- function(object, ...)
    stop("no 'ordicomparisons' method for objects of class ",
         paste(class(object), collapse = "/"))

#' Construct reference profiles for ordination contrasts
#'
#' Construct adjusted factor-level profiles from which ordination contrasts
#' are formed, optionally within groups. Each displayed row retains its
#' averaged model-matrix representation for use by `ordicomparisons()`.
#' The rows identify the quantities to combine in contrasts; they are not
#' individually zero-sum contrasts. Use `ordicomparisons()` directly for
#' routine comparisons and this constructor to inspect or name profile rows.
#'
#' @param object A fitted constrained ordination model.
#' @param variable The name of one focal factor in the model.
#' @param newdata The population over which profiles are averaged. `NULL`
#'   uses the analysis data, `"balanced"` uses a balanced reference grid,
#'   `"mean"` uses one mean-or-mode profile, and a data frame supplies an
#'   explicit population.
#' @param by Optional model variables defining groups within which profiles
#'   are averaged.
#' @param wts Optional finite, non-negative averaging weights or the name of
#'   a numeric weight column in `newdata`. Each group needs a positive weight.
#' @param ... Arguments passed to methods.
#'
#' @return A data frame of class `vegan_contrasts` with a `Profile` column
#'   identifying its rows. The `profile.matrix` attribute retains the averaged
#'   design rows; other attributes retain the model, grouping, reference
#'   population, weights, and extrapolation diagnostics.
#' @export
ordicontrasts <- function(object, ...)
    UseMethod("ordicontrasts")

#' Default comparison-profiles method
#'
#' Report that comparison profiles are unavailable for an object.
#'
#' @param object An unsupported object.
#' @param ... Unused arguments.
#'
#' @return This method always throws an error.
#' @keywords internal
ordicontrasts.default <- function(object, ...)
    stop("no 'ordicontrasts' method for objects of class ",
         paste(class(object), collapse = "/"))

#' Comparison profiles for constrained ordination
#'
#' Rebuild the fitted model matrix under every level of one factor and average
#' it over an observed, balanced, mean-or-mode, or user-supplied reference
#' population. This method is inherited by `rda`, `cca`, `dbrda`, and retained
#' `adonis2` fits.
#' Each profile represents a factor level with the other variables averaged
#' or held fixed. `ordicomparisons()` tests differences between their design
#' rows using the fitted model's multivariate geometry.
#'
#' @inheritParams ordicontrasts
#'
#' @return A `vegan_contrasts` data frame with stable `b1`, `b2`, ... row
#'   identifiers, focal and grouping values, reference-row counts, and
#'   extrapolation flags.
#' @export
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
    if (!variable %in% all.vars(.comparisonFormula(object)))
        stop(gettextf("variable '%s' is not a constraint in the model",
                      variable))
    if (!is.null(by)) {
        if (!is.character(by) || any(!by %in% names(dat)))
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
    ## The fitted QR includes conditions before constraints, even when
    ## no ordination axes or conditioning inertia were retained.
    design <- qr.X(Q, ncol = length(Q$pivot))
    if (!is.null(object$rowsum))
        design <- design / sqrt(object$rowsum)
    nz <- ncol(design) - p
    Z <- if (nz) design[, seq_len(nz), drop = FALSE] else NULL
    X <- design[, nz + seq_len(p), drop = FALSE]

    base <- .comparisonNewdata(object, dat, variable, by, newdata)
    meta <- .profileMetadata(object, dat, base, variable, by, newdata, wts)
    base$..vegan_rowid <- seq_len(nrow(base))
    split <- .comparisonSplit(base, by)
    lev <- levels(dat[[variable]])
    ordered <- isTRUE(object$terminfo$ordered[[variable]])
    refs <- lapply(seq_along(split), function(g) {
        gd <- split[[g]]
        gw <- .comparisonWeights(gd, wts, nrow(base))
        R <- .comparisonRows(object, gd, variable, lev, gw, colnames(X))
        group <- if (length(by)) gd[1L, by, drop = FALSE] else NULL
        extra <- .comparisonExtrapolation(gd, dat, variable)
        list(reference.rows = R, group = group, weights = gw,
             n.reference = nrow(gd), extrapolated = any(extra),
             extrapolation = extra)
    })

    tabs <- lapply(refs, function(z) {
        tab <- data.frame(Profile = rep.int("", length(lev)),
                          stringsAsFactors = FALSE)
        tab[[variable]] <- factor(lev, levels = lev, ordered = ordered)
        if (!is.null(z$group))
            tab <- cbind(tab, z$group[rep.int(1L, length(lev)), , drop = FALSE],
                         row.names = NULL)
        tab$N.reference <- rep.int(z$n.reference, length(lev))
        tab$Extrapolated <- rep.int(z$extrapolated, length(lev))
        tab
    })
    out <- do.call(rbind, tabs)
    rownames(out) <- NULL
    out$Profile <- paste0("b", seq_len(nrow(out)))
    R <- do.call(rbind, lapply(refs, `[[`, "reference.rows"))
    rownames(R) <- out$Profile
    group.index <- rep(seq_along(refs), each = length(lev))

    attr(out, "profile.matrix") <- R
    attr(out, "group.index") <- group.index
    attr(out, "reference") <- refs
    attr(out, "reference.data") <- base[, names(base) != "..vegan_rowid",
                                         drop = FALSE]
    attr(out, "model") <- object
    attr(out, "constraints") <- X
    attr(out, "conditions") <- Z
    attr(out, "variable") <- variable
    attr(out, "levels") <- lev
    attr(out, "by") <- by
    attr(out, "newdata") <- newdata
    attr(out, "weights") <- wts
    attr(out, "population") <- meta$population
    attr(out, "averaged") <- meta$averaged
    attr(out, "fixed") <- meta$fixed
    attr(out, "weight.label") <- meta$weights
    attr(out, "extrapolation") <- .comparisonExtrapolation(base, dat,
                                                             variable)
    class(out) <- c("vegan_contrasts", "data.frame")
    out
}

#' Print comparison profiles
#'
#' Print profile rows together with the reference population, averaging,
#' fixed-value, grouping, and weight choices that define them. Identifiers
#' `b1`, `b2`, ... can be used in `ordicomparisons()` expressions.
#'
#' @param x A `vegan_contrasts` object.
#' @param ... Arguments passed to `print.data.frame()`.
#'
#' @return `x`, invisibly.
#' @export
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

#' Comparisons for constrained ordination
#'
#' Construct comparison profiles and pass them to the profile comparison method.
#' This is the one-call interface inherited by `rda`, `cca`, `dbrda`, and
#' retained `adonis2` fits.
#'
#' @inheritParams ordicontrasts
#' @inheritParams ordicomparisons
#'
#' @return A data frame of class `vegan_comparisons` with the comparison
#'   degrees of freedom, inertia, pseudo-F statistic, and permutation p-value.
#' @export
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

#' Compare fitted profiles
#'
#' Translate recipes or linear expressions into differences between rows of a
#' `vegan_contrasts` table, then test their additional constrained inertia
#' using contrast-specific reduced-model permutation pseudo-F tests. Profile
#' identifiers such as `b2 - b1 = 0` refer directly to the displayed rows, not
#' to previously calculated comparisons or test statistics.
#'
#' @param object A `vegan_contrasts` object.
#' @inheritParams ordicomparisons
#'
#' @return A data frame of class `vegan_comparisons`.
#' @export
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
    perm <- getPermuteMatrix(permutations, nrow(basis), strata = strata)

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
            group <- if (length(by)) object[take[1L], by, drop = FALSE]
                     else NULL
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
        if (!is.null(spec$group)) {
            tab <- cbind(spec$group, tab, row.names = NULL)
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

#' Print constrained-ordination comparisons
#'
#' Print the comparison table with a short description of the fitted method.
#'
#' @param x A `vegan_comparisons` object.
#' @param ... Arguments passed to `print.data.frame()`.
#'
#' @return `x`, invisibly.
#' @export
print.vegan_comparisons <- function(x, ...)
{
    cat("Permutation comparisons for constrained ordination\n")
    cat("Variable:", attr(x, "variable"), "\n\n")
    print.data.frame(x, ...)
    invisible(x)
}

#' Recover analysis data for comparisons
#'
#' Evaluate the original data call and apply the fitted subset and missing-data
#' action. Restore the fitted factor levels after row selection. Raw data are
#' needed so transformed terms can be rebuilt faithfully.
#'
#' @param object A fitted constrained ordination model.
#'
#' @return A data frame aligned with the observations in `object`.
#' @keywords internal
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

#' Construct the comparison population
#'
#' Select the observed population, build a balanced or mean-or-mode reference
#' grid, or validate explicit new data.
#'
#' @param object A fitted constrained ordination model.
#' @param dat Analysis data.
#' @param variable Focal factor name.
#' @param by Optional grouping variables.
#' @param newdata Requested comparison population.
#'
#' @return A data frame used to construct counterfactual model matrices.
#' @keywords internal
.comparisonNewdata <- function(object, dat, variable, by, newdata)
{
    if (is.null(newdata))
        return(dat)
    if (is.character(newdata)) {
        choices <- c("balanced", "mean")
        if (length(newdata) != 1L || !newdata %in% choices)
            stop("character 'newdata' must be \"balanced\" or \"mean\"")
        vars <- unique(c(setdiff(all.vars(.comparisonFormula(object)),
                                 variable), by))
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

#' Select a typical value
#'
#' Return the arithmetic mean of a numeric variable or the mode of another
#' variable. Factor ties are resolved by fitted level order; other ties are
#' resolved by order of first appearance.
#'
#' @param x A model variable.
#'
#' @return A length-one value retaining the class and levels of `x` where
#'   applicable.
#' @keywords internal
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

#' Construct mean-or-mode reference profiles
#'
#' Construct one typical reference row globally, or one within every observed
#' combination of `by` variables. Numeric variables are represented by their
#' means and other non-grouping variables by their modes.
#'
#' @param dat Analysis data.
#' @param vars Model variables required in the reference data.
#' @param by Optional grouping variables.
#'
#' @return A data frame containing one mean-or-mode profile per group.
#' @keywords internal
.comparisonMeanData <- function(dat, vars, by)
{
    if (!length(vars))
        return(dat[1L, vars, drop = FALSE])
    groups <- if (length(by)) .comparisonSplit(dat, by) else list(dat)
    rows <- lapply(groups, function(gd) {
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

#' Describe profile averaging
#'
#' Summarize the reference population, averaged variables, fixed numeric or
#' user-supplied values, and weights used to construct comparison profiles. The
#' metadata is printed with the profile table and retained for diagnostics.
#'
#' @param object A fitted constrained ordination model.
#' @param dat Analysis data.
#' @param base Reference data used for profile construction.
#' @param variable Focal factor name.
#' @param by Optional grouping variables.
#' @param newdata Original `newdata` specification.
#' @param wts Original averaging-weight specification.
#'
#' @return A list of printable metadata.
#' @keywords internal
.profileMetadata <- function(object, dat, base, variable, by, newdata, wts)
{
    vars <- setdiff(all.vars(.comparisonFormula(object)), c(variable, by))
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

#' Diagnose extrapolation in comparison data
#'
#' Identify numeric values outside the range used to fit the model. Factor
#' levels are validated later by `model.matrix()` and its fitted `xlev`.
#'
#' @param newdata Comparison population.
#' @param dat Analysis data.
#' @param variable Focal factor name.
#'
#' @return A named logical vector for numeric model variables.
#' @keywords internal
.comparisonExtrapolation <- function(newdata, dat, variable)
{
    vars <- intersect(setdiff(names(newdata), variable), names(dat))
    vars <- vars[vapply(dat[vars], is.numeric, logical(1))]
    setNames(vapply(vars, function(v) {
        r <- range(dat[[v]], na.rm = TRUE)
        any(newdata[[v]] < r[1L] | newdata[[v]] > r[2L], na.rm = TRUE)
    }, logical(1)), vars)
}

#' Split a comparison population by moderators
#'
#' @param dat Comparison population.
#' @param by Optional grouping variables.
#'
#' @return A named list of data frames.
#' @keywords internal
.comparisonSplit <- function(dat, by)
{
    if (!length(by))
        return(list(dat))
    split(dat, interaction(dat[by], drop = TRUE, lex.order = TRUE))
}

#' Resolve comparison averaging weights
#'
#' @param dat Comparison population.
#' @param wts `NULL`, a numeric vector, or a weight-column name.
#' @param n Total rows before optional `by` splitting.
#'
#' @return A normalised numeric vector.
#' @keywords internal
.comparisonWeights <- function(dat, wts, n)
{
    if (is.null(wts))
        w <- rep.int(1, nrow(dat))
    else if (is.character(wts) && length(wts) == 1L && wts %in% names(dat))
        w <- dat[[wts]]
    else if (is.numeric(wts) && length(wts) == n)
        w <- wts[dat$..vegan_rowid]
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

#' Construct the constraint formula
#'
#' Remove `Condition()` terms while preserving ordinary main effects,
#' interactions, and transformations.
#'
#' @param object A fitted constrained ordination model.
#'
#' @return A one-sided formula for the constraint model matrix.
#' @keywords internal
.comparisonFormula <- function(object)
{
    Terms <- terms(object)
    labs <- attr(Terms, "term.labels")
    labs <- labs[!grepl("^Condition\\(", labs)]
    reformulate(labs, intercept = attr(Terms, "intercept"),
                env = environment(Terms))
}

#' Construct averaged model-matrix rows
#'
#' Set the focal factor to each fitted level, reconstruct the complete model
#' matrix, and average over the selected population. The fitted terms retain
#' prediction calls supplied by `makepredictcall()`, including the original
#' polynomial coefficients, spline knots and scaling constants.
#'
#' @param object A fitted constrained ordination model.
#' @param dat Comparison population.
#' @param variable Focal factor name.
#' @param levels Fitted factor levels.
#' @param w Normalized averaging weights.
#' @param columns Constraint-matrix column names.
#'
#' @return A levels-by-coefficients reference-row matrix.
#' @keywords internal
.comparisonRows <- function(object, dat, variable, levels, w, columns)
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
    xlev <- xlev[names(xlev) %in% all.vars(Terms)]
    ordered <- isTRUE(object$terminfo$ordered[[variable]])
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

#' Detect comparisons written against profile rows
#'
#' Distinguish global expressions or numeric weights over displayed `b` rows
#' from recipes and factor-level templates that should be repeated within each
#' `by` group.
#'
#' @param comparison Comparison specification.
#' @param ids Profile identifiers.
#' @param nlevels Number of focal-factor levels per group.
#'
#' @return A single logical value.
#' @keywords internal
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

#' Generate comparison weights
#'
#' Translate built-in recipes, equations, lists, vectors, and matrices into
#' weights on factor-level reference rows.
#'
#' @param comparison Comparison specification.
#' @param keys Names of the quantities being compared.
#' @param reference Optional reference quantity.
#'
#' @return A named contrast matrix.
#' @keywords internal
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
            if (is.na(ref)) stop("'reference' is not one of the quantities")
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
        A <- do.call(rbind, lapply(comparison, function(z)
            .comparisonContrasts(z, keys, reference)))
        if (!is.null(names(comparison)) && all(nzchar(names(comparison))))
            rownames(A) <- names(comparison)
    } else if (is.numeric(comparison)) {
        if (is.matrix(comparison)) {
            A <- comparison
        } else {
            A <- matrix(comparison, nrow = 1L,
                        dimnames = list(NULL, names(comparison)))
        }
        if (!is.null(colnames(A))) {
            if (!all(keys %in% colnames(A)))
                stop("contrast matrix columns must contain all quantities")
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

#' Compile differences between fitted profiles
#'
#' Normalise each row of weights and apply it to differences from a
#' common profile. A shared offset in the profiles then cancels before
#' multiplication. Changing the overall weight scale leaves the null
#' hypothesis unchanged.
#'
#' @param A Comparison weights with one column per profile.
#' @param R Matrix of averaged model-matrix rows.
#'
#' @return A matrix with one row per comparison.
#' @seealso .comparisonContrasts, .comparisonTest
#' @keywords internal
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

#' Parse a safe linear comparison expression
#'
#' @param expression A character equation.
#' @param keys Names allowed in the equation.
#'
#' @return A numeric vector of linear weights.
#' @keywords internal
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

#' Label comparison rows
#'
#' @param A Comparison matrix.
#' @param keys Quantity names.
#'
#' @return Character labels describing non-zero weights.
#' @keywords internal
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

#' Test one compiled comparison
#'
#' Separate the null and tested spaces in the fitted design's orthonormal
#' coordinates. Use QR projections to calculate inertia before passing the
#' test to the reduced-model permutation machinery.
#'
#' @param object Original constrained ordination model.
#' @param basis Orthonormal basis of the fitted design, with conditions first.
#' @param triangular Constraint block of the fitted triangular QR factor.
#' @param pivot Order of the constraint columns in `triangular`.
#' @param L Compiled comparison matrix.
#' @param permutations Shared permutation matrix.
#' @param parallel Parallel-processing specification.
#'
#' @return A list of test statistics.
#' @seealso .comparisonMatrix, permutest.cca
#' @keywords internal
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

#' Convert comparisons to an ANOVA table
#'
#' @param x A `vegan_comparisons` result.
#' @param object Original constrained ordination model.
#'
#' @return An `anova.cca` data frame.
#' @keywords internal
.comparisonsAnova <- function(x, object)
{
    keep <- x[, c("Df", "Inertia", "F", "Pr(>F)"), drop = FALSE]
    varname <- if (inherits(object, c("capscale", "dbrda")) &&
                   object$adjust == 1)
        "SumOfSqs" else if (inherits(object, "rda")) "Variance" else "ChiSquare"
    names(keep)[2L] <- varname
    rownames(keep) <- make.unique(x$Comparison)
    attr(keep, "heading") <- c("Permutation test for specified comparisons\n",
                               paste("Model:", c(object$call)))
    attr(keep, "Random.seed") <- attr(x, "Random.seed")
    attr(keep, "control") <- attr(x, "control")
    class(keep) <- c("anova.cca", "anova", "data.frame")
    keep
}
