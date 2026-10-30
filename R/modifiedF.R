## Components of the modified pseudo-F statistic of Anderson et al.
## (2017).  These are kept separate from the permutation machinery so
## that other error terms and variance-component coefficients can be
## supplied here if constrained ordination gains random-factor models.

`modifiedFGroups` <-
    function(x, groups)
{
    if (inherits(x, "capscale"))
        stop("'permustat = \"modified-F\"' is not available for capscale")
    if (!inherits(x, "rda"))
        stop("'permustat = \"modified-F\"' is not available for cca; ",
             "use adonis2() or dbrda() with explicit dissimilarities")
    if (!is.null(x$pCCA) && x$pCCA$rank > 0)
        stop("'permustat = \"modified-F\"' is not available with Conditions")
    if (is.null(groups))
        stop("'groups' must identify the groups with heterogeneous dispersions")

    ## Exact equality of projection spaces is necessary but not sufficient
    ## for the currently supported design: a numeric covariate can happen to
    ## encode the same space. Limit this implementation to declared fixed
    ## factors, leaving quantitative and mixed/random terms to their own
    ## denominator and permutation designs.
    terminfo <- x$terminfo
    term.variables <- if (!is.null(terminfo))
        rownames(attr(terminfo$terms, "factors"))
    else character(0)
    if (!length(term.variables) ||
        any(!term.variables %in% names(terminfo$xlev)))
        stop("all constraint terms must be factors for the modified-F test")

    n <- nobs(x)
    if (length(groups) != n)
        stop(gettextf("'groups' has length %d, but the model has %d observations",
                      length(groups), n))
    if (anyNA(groups))
        stop("missing values are not allowed in 'groups'")
    groups <- droplevels(as.factor(groups), exclude = NA)
    ngroup <- nlevels(groups)
    if (ngroup < 2L)
        stop("'groups' must define at least two groups")
    group.size <- as.vector(table(groups))
    names(group.size) <- levels(groups)
    if (any(group.size < 2L))
        stop("each level of 'groups' must contain at least two observations")

    ## The published statistic is an omnibus one-way test.  A saturated
    ## factorial model is admissible only when its constraint space is
    ## exactly the space spanned by the combined cell factor.
    Xgroup <- model.matrix(~ groups)[, -1L, drop = FALSE]
    Xgroup <- scale(Xgroup, scale = FALSE)
    Qgroup <- qr(Xgroup)
    Qmodel <- x$CCA$QR
    qrank <- x$CCA$qrank
    if (qrank != ngroup - 1L || Qgroup$rank != ngroup - 1L)
        stop("the model must be an omnibus one-way or saturated cell model for 'groups'")

    Xmodel <- qr.X(Qmodel, ncol = length(Qmodel$pivot))
    in.span <- function(Q, X) {
        z <- qr.resid(Q, X)
        sqrt(sum(z^2)) <= 1e-7 * max(1, sqrt(sum(X^2)))
    }
    if (!in.span(Qmodel, Xgroup) || !in.span(Qgroup, Xmodel))
        stop("the model constraints and the centred 'groups' indicators ",
             "must span the same space")

    list(groups = groups, size = group.size,
         index = split(seq_len(n), groups), rank = qrank)
}

## Build an error term from variance components. Keeping this interface
## independent of the response projection leaves room for future
## expected-mean-square implementations to provide different components,
## coefficients, degrees of freedom, and permutation units.
`modifiedFError` <-
    function(group.inertia, group.size,
             coefficients = 1 - group.size / sum(group.size))
{
    dispersion <- group.inertia / (group.size - 1)
    tol <- sqrt(.Machine$double.eps) * max(1, sum(abs(group.inertia)))
    if (any(dispersion < -tol))
        stop("negative within-group dispersion in modified pseudo-F calculation")
    dispersion[dispersion < 0] <- 0

    components <- coefficients * dispersion
    denominator <- sum(components)
    den.df <- if (denominator > 0) {
        denominator^2 / sum(components^2 / (group.size - 1))
    } else {
        NA_real_
    }
    list(denominator = denominator, df = den.df,
         dispersion = dispersion, components = components,
         coefficients = coefficients)
}

`modifiedFStatistic` <-
    function(E, Q, group.info, isDB)
{
    if (isDB) {
        fitted <- qr.fitted(Q, E)
        fitted <- qr.fitted(Q, t(fitted))
        residual <- qr.resid(Q, E)
        residual <- qr.resid(Q, t(residual))
        numerator <- sum(diag(fitted))
        row.inertia <- diag(residual)
    } else {
        fitted <- qr.fitted(Q, E)
        residual <- qr.resid(Q, E)
        numerator <- sum(fitted^2)
        row.inertia <- rowSums(residual^2)
    }

    group.inertia <- vapply(group.info$index,
                            function(i) sum(row.inertia[i]), numeric(1))
    error <- modifiedFError(group.inertia, group.info$size)
    statistic <- numerator / error$denominator

    list(numerator = numerator, denominator = error$denominator,
         statistic = statistic, dispersion = error$dispersion,
         coefficients = error$coefficients, components = error$components,
         den.df = error$df,
         residual = sum(group.inertia))
}

`modifiedFPermutations` <-
    function(indx, E, Q, group.info, isDB)
{
    if (!is.matrix(indx))
        indx <- matrix(indx, nrow = 1L)
    out <- matrix(NA_real_, nrow(indx), 3L,
                  dimnames = list(NULL, c("num", "den", "F")))
    for (i in seq_len(nrow(indx))) {
        take <- indx[i, ]
        Eperm <- if (isDB) E[take, take, drop = FALSE]
                 else E[take, , drop = FALSE]
        stat <- modifiedFStatistic(Eperm, Q, group.info, isDB)
        out[i, ] <- c(stat$numerator, stat$denominator, stat$statistic)
    }
    out
}
