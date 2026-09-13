#' Comparisons for constrained ordination models
#'
#' Construct scientifically meaningful comparisons of factor levels and test
#' them with contrast-specific permutation pseudo-F statistics.
#'
#' @param object A fitted constrained ordination model.
#' @param variables The name of one factor in the model.
#' @param comparison A built-in comparison name, character equations, or a
#'   numeric contrast vector or matrix.
#' @param reference Optional reference or control level. The first fitted
#'   factor level is used when this is omitted.
#' @param newdata The population over which comparisons are averaged. `NULL`
#'   uses the analysis data, `"balanced"` uses a balanced reference grid, and
#'   a data frame supplies an explicit population.
#' @param by Optional variables defining groups within which comparisons are
#'   averaged.
#' @param hypothesis Optional linear hypotheses combining the generated
#'   comparisons.
#' @param joint Test all resulting comparisons jointly.
#' @param permutations A permutation count, matrix, or `permute::how()` object.
#' @param strata Optional permutation strata.
#' @param parallel Number of parallel workers or a socket cluster.
#' @param adjust A method accepted by `stats::p.adjust()`.
#' @param wts Optional averaging weights or the name of a weight column.
#' @param ... Arguments passed to methods.
#'
#' @return A data frame of class `vegan_comparisons`.
#' @export
comparisons <- function(object, ...)
    UseMethod("comparisons")

#' Default comparisons method
#'
#' Report that no comparisons method is available for an object.
#'
#' @param object An unsupported object.
#' @param ... Unused arguments.
#'
#' @return This method always throws an error.
#' @keywords internal
comparisons.default <- function(object, ...)
    stop("no 'comparisons' method for objects of class ",
         paste(class(object), collapse = "/"))

#' Comparisons for constrained ordination
#'
#' Compile factor-level comparisons into model-matrix hypotheses and test
#' their additional constrained inertia under a contrast-specific reduced
#' permutation model. This method is inherited by `rda`, `cca`, and `dbrda`
#' objects.
#'
#' @inheritParams comparisons
#'
#' @return A data frame of class `vegan_comparisons` with the hypothesis
#'   degrees of freedom, inertia, pseudo-F statistic, and permutation p-value.
#' @export
comparisons.cca <-
    function(object, variables, comparison = "pairwise", reference = NULL,
             newdata = NULL, by = NULL, hypothesis = NULL, joint = FALSE,
             permutations = how(nperm = 999), strata = NULL,
             parallel = getOption("mc.cores"), adjust = "none", wts = NULL,
             ...)
{
    if (is.null(object$terms) || is.null(object$terminfo))
        stop("comparisons need a model fitted with the formula interface")
    if (length(variables) != 1L || !is.character(variables))
        stop("'variables' must name exactly one factor")
    variable <- variables[[1L]]
    dat <- .comparisonData(object)
    if (!variable %in% names(dat))
        stop(gettextf("variable '%s' was not found in the model data",
                      variable))
    if (!is.factor(dat[[variable]]))
        stop(gettextf("variable '%s' must be a factor", variable))
    if (!is.null(by)) {
        if (!is.character(by) || any(!by %in% names(dat)))
            stop("'by' must name variables in the model data")
        if (variable %in% by)
            stop("the focal variable cannot also be a 'by' variable")
    }
    adjust <- match.arg(adjust, p.adjust.methods)
    if (length(joint) != 1L || is.na(joint))
        stop("'joint' must be TRUE or FALSE")

    Xs <- model.matrix(object)
    if (is.list(Xs)) {
        Z <- Xs$Conditions
        X <- Xs$Constraints
    } else {
        Z <- NULL
        X <- Xs
    }
    if (is.null(X) || !ncol(X))
        stop("model has no constraints")
    if (object$CCA$qrank < ncol(X))
        stop("comparisons are not available for aliased constraint matrices")

    base <- .comparisonNewdata(object, dat, variable, newdata)
    extrapolation <- .comparisonExtrapolation(base, dat, variable)
    base$..vegan_rowid <- seq_len(nrow(base))
    split <- .comparisonSplit(base, by)
    lev <- levels(dat[[variable]])
    perm <- getPermuteMatrix(permutations, nobs(object), strata = strata)
    refs <- lapply(seq_along(split), function(g) {
        gd <- split[[g]]
        gw <- .comparisonWeights(gd, wts, nrow(base))
        R <- .comparisonRows(object, gd, variable, lev, gw, colnames(X))
        A <- .comparisonContrasts(comparison, lev, reference)
        group <- if (length(by)) gd[1L, by, drop = FALSE] else NULL
        list(reference.rows = R, contrast.weights = A, group = group)
    })

    ## A hypothesis with 'by' combines the complete ordered set of grouped
    ## comparisons. This permits interaction contrasts such as differences
    ## of simple effects. Without a hypothesis, groups remain separate rows.
    if (length(by) && !is.null(hypothesis)) {
        L0 <- do.call(rbind, lapply(refs, function(z) {
            a <- z$contrast.weights %*% z$reference.rows
            grp <- paste(paste(by, unlist(z$group), sep = "="),
                         collapse = ", ")
            rownames(a) <- paste0(rownames(a), " | ", grp)
            a
        }))
        ident <- diag(nrow(L0))
        rownames(ident) <- colnames(ident) <- rownames(L0)
        H <- .comparisonHypothesis(ident, hypothesis)
        L <- H %*% L0
        specs <- lapply(seq_len(nrow(L)), function(i)
            list(label = rownames(L)[i], L = L[i, , drop = FALSE],
                 group = NULL))
    } else {
        specs <- unlist(lapply(refs, function(z) {
            A <- z$contrast.weights
            if (!is.null(hypothesis))
                A <- .comparisonHypothesis(A, hypothesis)
            if (joint && nrow(A) > 1L) {
                return(list(list(label = paste(rownames(A), collapse = "; "),
                                      L = A %*% z$reference.rows,
                                      group = z$group)))
            }
            lapply(seq_len(nrow(A)), function(i)
                list(label = rownames(A)[i],
                     L = A[i, , drop = FALSE] %*% z$reference.rows,
                     group = z$group))
        }), recursive = FALSE)
    }
    if (joint && length(specs) > 1L) {
        L <- do.call(rbind, lapply(specs, `[[`, "L"))
        specs <- list(list(label = paste(vapply(specs, `[[`, "", "label"),
                                             collapse = "; "),
                           L = L, group = NULL))
    }

    ans <- lapply(specs, function(spec) {
        z <- .comparisonTest(object, X, Z, spec$L, perm, parallel)
        tab <- data.frame(Comparison = spec$label, Df = z$df,
                          Inertia = z$inertia,
                          R2 = z$inertia / object$tot.chi,
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
    attr(out, "comparison.matrices") <- list(reference = refs,
                                               hypotheses = lapply(specs, `[[`,
                                                                   "L"))
    attr(out, "model") <- object
    attr(out, "estimable") <- rep.int(TRUE, nrow(out))
    attr(out, "extrapolation") <- extrapolation
    attr(out, "variable") <- variable
    attr(out, "comparison") <- comparison
    attr(out, "hypothesis") <- hypothesis
    attr(out, "newdata") <- newdata
    attr(out, "by") <- by
    attr(out, "adjust") <- adjust
    attr(out, "weights") <- wts
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
#' action. Raw data are needed so transformed terms can be rebuilt faithfully.
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
    dat <- eval(call$data, environment(call$formula), .GlobalEnv)
    if (inherits(dat, "environment"))
        dat <- as.data.frame(as.list(dat))
    else
        dat <- as.data.frame(dat)
    if (!is.null(object$subset))
        dat <- dat[object$subset, , drop = FALSE]
    if (!is.null(object$na.action))
        dat <- dat[-as.integer(object$na.action), , drop = FALSE]
    if (nrow(dat) != nobs(object))
        stop("could not align reconstructed model data with fitted ",
             "observations")
    dat
}

#' Construct the comparison population
#'
#' Select the observed population, build a balanced reference grid, or validate
#' explicit new data.
#'
#' @param object A fitted constrained ordination model.
#' @param dat Analysis data.
#' @param variable Focal factor name.
#' @param newdata Requested comparison population.
#'
#' @return A data frame used to construct counterfactual model matrices.
#' @keywords internal
.comparisonNewdata <- function(object, dat, variable, newdata)
{
    if (is.null(newdata))
        return(dat)
    if (is.character(newdata)) {
        if (length(newdata) != 1L || newdata != "balanced")
            stop("character 'newdata' must be \"balanced\"")
        vars <- setdiff(all.vars(.comparisonFormula(object)), variable)
        vals <- lapply(vars, function(v) {
            x <- dat[[v]]
            if (is.factor(x)) levels(x)
            else if (is.numeric(x)) mean(x)
            else unique(x)
        })
        names(vals) <- vars
        if (!length(vals))
            return(dat[1L, , drop = FALSE])
        return(expand.grid(vals, KEEP.OUT.ATTRS = FALSE,
                           stringsAsFactors = FALSE))
    }
    if (!is.data.frame(newdata))
        stop("'newdata' must be NULL, \"balanced\", or a data frame")
    newdata
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
#' @return A normalized numeric vector.
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
    if (anyNA(w) || any(w < 0) || !sum(w) > 0)
        stop("comparison weights must be non-negative, finite, and sum ",
             "above zero")
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
    labs <- attr(terms(object), "term.labels")
    labs <- labs[!grepl("^Condition\\(", labs)]
    reformulate(labs, intercept = attr(terms(object), "intercept"),
                env = environment(formula(object)))
}

#' Construct averaged model-matrix rows
#'
#' Set the focal factor to each fitted level, reconstruct the complete model
#' matrix, and average over the selected population.
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
    f <- .comparisonFormula(object)
    ctr <- object$terminfo$contrasts
    xlev <- object$terminfo$xlev
    xlev <- xlev[names(xlev) %in% all.vars(f)]
    ordered <- isTRUE(object$terminfo$ordered[[variable]])
    rows <- lapply(levels, function(lev) {
        nd <- dat
        nd[[variable]] <- factor(lev, levels = levels,
                                 ordered = ordered)
        mm <- model.matrix(f, nd, contrasts.arg = ctr,
                           xlev = xlev)
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

#' Generate comparison weights
#'
#' Translate built-in recipes, equations, lists, vectors, and matrices into
#' weights on factor-level reference rows.
#'
#' @param comparison Comparison specification.
#' @param keys Names of the quantities being compared.
#' @param reference Optional reference quantity.
#' @param zero Require weights to sum to zero.
#'
#' @return A named contrast matrix.
#' @keywords internal
.comparisonContrasts <- function(comparison, keys, reference = NULL,
                                 zero = TRUE)
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
            if (is.na(ref)) stop("'reference' is not a fitted factor level")
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
            .comparisonContrasts(z, keys, reference, zero)))
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
                stop("contrast matrix columns must contain all factor levels")
            A <- A[, keys, drop = FALSE]
        } else if (ncol(A) != k) {
            stop("numeric comparisons must have one column per factor level")
        }
        if (is.null(rownames(A))) rownames(A) <- .comparisonLabels(A, keys)
    } else {
        stop("unsupported 'comparison' specification")
    }
    colnames(A) <- keys
    if (zero && any(abs(rowSums(A)) > sqrt(.Machine$double.eps)))
        stop("comparison weights must sum to zero")
    A
}

#' Compose generated comparisons into new hypotheses
#'
#' @param A First-stage comparison matrix.
#' @param hypothesis Character equations or a numeric matrix.
#'
#' @return A composed comparison matrix.
#' @keywords internal
.comparisonHypothesis <- function(A, hypothesis)
{
    keys <- paste0("b", seq_len(nrow(A)))
    recipe <- length(hypothesis) == 1L && is.character(hypothesis) &&
        hypothesis %in% c("pairwise", "sequential", "consecutive",
                          "reference", "trt.vs.ctrl", "deviation",
                          "effect", "helmert", "poly")
    if (is.character(hypothesis) && !recipe) {
        labs <- rownames(A)
        if (anyDuplicated(labs))
            labs <- rep.int("", length(labs))
        allowed <- c(keys, setdiff(labs[nzchar(labs)], keys))
        C <- do.call(rbind, lapply(hypothesis, .comparisonExpression,
                                   keys = allowed))
        colnames(C) <- allowed
        B <- C[, keys, drop = FALSE]
        for (i in seq_along(labs))
            if (nzchar(labs[i]) && labs[i] %in% colnames(C))
                B[, i] <- B[, i] + C[, labs[i]]
        rownames(B) <- if (!is.null(names(hypothesis)) &&
                           all(nzchar(names(hypothesis)))) names(hypothesis)
                      else hypothesis
    } else {
        B <- .comparisonContrasts(hypothesis, keys, zero = FALSE)
    }
    out <- B %*% A
    rownames(out) <- rownames(B)
    out
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
    if (abs(z$constant) > sqrt(.Machine$double.eps))
        stop("only hypotheses with a zero right-hand side are supported")
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
#' Construct null and tested bases, fit a temporary partial constrained model,
#' and obtain its reduced-model permutation pseudo-F statistic.
#'
#' @param object Original constrained ordination model.
#' @param X Constraint matrix.
#' @param Z Original condition matrix.
#' @param L Compiled hypothesis matrix.
#' @param permutations Shared permutation matrix.
#' @param parallel Parallel-processing specification.
#'
#' @return A list of test statistics.
#' @keywords internal
.comparisonTest <- function(object, X, Z, L, permutations, parallel)
{
    tol <- sqrt(.Machine$double.eps)
    if (!nrow(L) || max(abs(L)) <= tol)
        stop("comparison has rank zero")
    ql <- qr(t(L), tol = tol)
    rank <- ql$rank
    Q <- qr.Q(ql, complete = TRUE)
    M <- Q[, seq_len(rank), drop = FALSE]
    N <- if (rank < ncol(X)) Q[, (rank + 1L):ncol(X), drop = FALSE]
         else matrix(numeric(), nrow = ncol(X), ncol = 0L)
    X1 <- X %*% M
    X0 <- if (ncol(N)) X %*% N else NULL
    Z0 <- if (is.null(Z)) X0 else if (is.null(X0)) Z else cbind(Z, X0)
    fit <- suppressMessages(ordConstrained(ordiYbar(object, "initial"),
                                           X = X1, Z = Z0, method = "pass"))
    if (is.null(fit$CCA) || fit$CCA$qrank != rank)
        stop("comparison is not estimable in the fitted model")
    if (is.null(fit$CA) || fit$CA$rank == 0)
        stop("comparison test has no residual degrees of freedom")
    fit$call <- object$call
    fit$adjust <- object$adjust
    tst <- permutest.cca(fit, permutations = permutations,
                         model = "reduced", parallel = parallel)
    eps <- sqrt(.Machine$double.eps)
    p <- (sum(tst$F.perm >= tst$F.0 - eps) + 1) / (tst$nperm + 1)
    list(df = rank, inertia = fit$CCA$tot.chi,
         F = unname(tst$F.0), p = p)
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
