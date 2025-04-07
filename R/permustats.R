### Functions to extract permutation statististic or null model
### results from various vegan objects.

## extract items as 'statistic' and 'permutations'. Specific methods
## towards the end of this file

`permustats` <-
    function(x, ...)
{
    UseMethod("permustats")
}

## something like str()
`print.permustats` <-
    function(x, ...)
{
    print(str(x))
    invisible(x)
}

### modelled after print.oecosimu (should perhaps have oecosimu() args
### like 'alternative'

`summary.permustats` <-
    function(object, interval = 0.95, alternative, ...)
{
    TAB <- c("two.sided", "greater", "less")
    if (missing(alternative))
        alt <- TAB[match(object$alternative, TAB)]
    else
        alt <- match.arg(alternative, TAB, several.ok = TRUE)
    if (any(is.na(alt)))
        stop("alternative missing")
    nstat <- length(object$statistic)
    nalt <- length(alt)
    ## Replicate alternative to length of statistic
    if ((nalt < nstat) && identical(nalt, 1L)) {
        alt <- rep(alt, length.out = nstat)
    }
    compint <- (1 - interval) / 2
    PROBS <- list(two.sided = c(compint, 0.5, interval + compint),
                  greater = c(NA, 0.5, interval),
                  less = c(1 - interval, 0.5, NA))
    probs <- PROBS[alt]
    ## take care that permutations are in a column matrix
    permutations <- as.matrix(object$permutations)
    object$means <- colMeans(permutations)
    sd <- apply(permutations, 2, sd)
    object$z <-
        (object$statistic - object$means)/sd
    qFun <- function(i, sim, probs) {
        quantile(sim[, i], probs = probs[[i]], na.rm = TRUE)
    }
    object$quantile <- lapply(seq_along(probs), qFun, sim = permutations, probs = probs)
    object$quantile <- do.call("rbind", object$quantile)
    dimnames(object$quantile) <- list(NULL, c("lower", "median", "upper"))
    object$interval <- interval
    ## P-values
    if (is.integer(object$statistic) && is.integer(permutations)) {
        pless <- rowSums(object$statistic >= t(permutations), na.rm = TRUE)
        pmore <- rowSums(object$statistic <= t(permutations), na.rm = TRUE)
    } else {
        EPS <- sqrt(.Machine$double.eps)
        pless <- rowSums(object$statistic + EPS >= t(permutations),
                         na.rm = TRUE)
        pmore <- rowSums(object$statistic - EPS <= t(permutations),
                         na.rm = TRUE)
    }
    nsimul <- nrow(permutations)
    if (any(is.na(permutations))) {
        warning("some simulated values were NA and were removed")
        nsimul <- nsimul - colSums(is.na(permutations))
    }
    p <- rep(NA, length(object$statistic))
    for(i in seq_along(p))
        p[i] <- switch(alt[i],
                       two.sided = 2*pmin.int(pless[i], pmore[i]),
                       greater = pmore[i],
                       less = pless[i])
    object$p <- pmin.int(1, (p + 1)/(nsimul + 1))
    ## out
    class(object) <- "summary.permustats"
    object
}

`print.summary.permustats` <- function(x, ...) {
    m <- cbind("statistic" = x$statistic,
               "SES" = x$z,
               "mean" = x$means,
               x$quantile,
               "Pr(perm)" = x$p)
    cat("\n")
    printCoefmat(m, tst.ind = 1:(ncol(m)-1), na.print = "", ...)
    writeLines(strwrap(paste0("(Interval (Upper - Lower) = ", x$interval, ")", sep = ""),
                       initial = "\n"))
    invisible(x)
}

### combine permustats objects. Function checks that statistic field
### is equal (name, value) before proceeding, sees if the alternative
### is equal, and then combines permutations.

`c.permustats` <-
    function(..., recursive = FALSE)
{
    mods <- list(...)
    ## check stats
    stats <- lapply(mods, function(z) z$statistic)
    if (!all(sapply(stats[-1], function(z) identical(stats[[1]], z))))
        stop("statistics are not equal")
    stats <- stats[[1]]
    ## check alternative
    alt <- lapply(mods, function(z) z$alternative)
    if (all(sapply(alt[-1], function(z) identical(alt[[1]], z))))
        alt <- alt[[1]]
    else
        alt <- NA
    ## combine permutations
    p <- do.call(rbind, lapply(mods, function(z) z$permutations))
    ## return permustats
    structure(list(statistic = stats,
                   permutations = p,
                   alternative = alt),
              class = "permustats")
}

### Wrapper function to allow calling S3 method functions densityplot
### & qqmath without loading (attaching) lattice package.

`permulattice` <-
    function(x, plot = c("densityplot", "qqmath"), observed = TRUE,
             axislab = "Permutations", ...)
{
    plot <- match.arg(plot)
    switch(plot,
           "densityplot" =
               densityplot(x, observed = observed, xlab = axislab, ...),
           "qqmath" =
               qqmath(x, observed = observed, ylab = axislab, ...)
           )
}

### lattice::densityplot

`densityplot.permustats` <-
    function(x, data, observed = TRUE, xlab = "Permutations", ...)
{
    obs <- x$statistic
    sim <- as.matrix(x$permutations)
    if (observed)
        sim <- rbind(x$statistic, sim)
    nm <- names(obs)[col(sim)]
    densityplot( ~ as.vector(sim) | factor(nm, levels = unique(nm)),
                xlab = xlab,
                panel = function(x, ...) {
                    panel.densityplot(x, ...)
                    if (observed)
                        panel.abline(v = obs[panel.number()], ...)
                },
                ...)
}

### simple density: normally densityplot should be used (or I suggest
### so), but we also offer basic density. This can be either with or
### without observed statistic.

`density.permustats` <-
    function(x, observed = TRUE, ...)
{
    ## only works with statistic
    if (length(x$statistic) > 1)
        stop(gettextf("only works with one statistic: you got %d",
                      length(x$statistic)))
    p <- x$permutations
    if (observed)
        p <- c(x$statistic, p)
    out <- density(p)
    out$call <- match.call()
    out$call[[1]] <- as.name("density")
    out
}

### QQ-plot against Guaussian distribution

`qqnorm.permustats` <-
    function(y, observed = TRUE, ...)
{
    ## only works with statistic
    if (length(y$statistic) > 1)
        stop(gettextf("only works with one statistic: you got %d",
                      length(y$statistic)))
    p <- y$permutations
    if (observed)
        p <- c(y$statistic, p)
    q <- qqnorm(p, ...)
    if (observed)
        abline(h = y$statistic, ...)
    invisible(q)
}

## lattice::qqmath

`qqmath.permustats` <-
    function(x, data, observed = TRUE, sd.scale = FALSE,
             ylab = "Permutations", ...)
{
    ## sd.scale: standardize before use
    if (sd.scale) {
        x$permutations <- scale(x$permutations)
        x$statistic <- (x$statistic - attr(x$permutations, "scaled:center"))/
            attr(x$permutations, "scaled:scale")
    }
    obs <- x$statistic
    if (observed)
        sim <- rbind(x$statistic, as.matrix(x$permutations))
    else
        sim <- as.matrix(x$permutations)
    nm <- names(obs)[col(sim)]
    qqmath( ~ as.vector(sim) | factor(nm, levels = unique(nm)),
                ylab = ylab,
                panel = function(x, ...) {
                    panel.qqmath(x, ...)
                    if (observed)
                        panel.abline(h = obs[panel.number()], ...)
                },
                ...)
}

## boxplot for (standardized) effect size: permutations are centred to
## the value of statistic, and optionally standardized to equal sd
## w.r.t. to column mean (not to the statistics). This shows the
## effect size, or the deviation from the observed statistic (= 0).

`boxplot.permustats` <-
    function(x, scale = FALSE, names, ...)
{
    p <- x$permutations
    if (isTRUE(scale))
        scale <- apply(p, 2, sd, na.rm = TRUE)
    p <- scale(p, center = x$statistic, scale = scale)
    if (missing(names))
        names <- attr(x$statistic, "names")
    boxplot(p, names = names, ...)
}

## pairs plot permuted variables against each other

`pairs.permustats` <-
    function(x, ...)
{
    p <- x$permutations
    colnames(p) <- attr(x$statistic, "names")
    pairs(p, ...)
}

###
### specific methods to extract permustats
###

`permustats.anosim` <-
    function(x, ...)
{
    structure(list(
        "statistic" = structure(x$statistic, names = "R"),
        "permutations" = x$perm,
        "alternative" = "greater"),
              class = "permustats")
}

## summary: see file R/summary.anosim.R

`permustats.mantel` <-
    function(x, ...)
{
    structure(list(
        "statistic" = structure(x$statistic, names = "r"),
        "permutations" = x$perm,
        "alternative" = "greater"),
              class = "permustats")
}

`summary.mantel` <-
    function(object, ...)
{
    summary(permustats(object, ...))
}

`permustats.mrpp` <-
    function(x, ...)
{
    structure(list(
        "statistic" = structure(x$delta, names = "delta"),
        "permutations" = x$boot.deltas,
        "alternative" = "less"),
              class = "permustats")
}

`summary.mrpp` <-
    function(object, ...)
{
    summary(permustats(object, ...))
}

`permustats.oecosimu` <-
    function(x, ...)
{
    structure(list(
        "statistic" = x$oecosimu$statistic,
        "permutations" = t(x$oecosimu$simulated),
        "alternative" = x$oecosimu$alternative),
              class = "permustats")
}

`summary.oecosimu` <-
    function(object, ...)
{
    summary(permustats(object, ...))
}

`permustats.ordiareatest` <-
    function(x, ...)
{
    structure(list(
        "statistic" = x$areas,
        "permutations" = t(x$permutations),
        "alternative" = "less"),
              class = "permustats")
}

`summary.ordiareatest` <-
    function(object, ...)
{
    summary(permustats(object, ...))
}

`permustats.permutest.cca` <-
    function(x, ...)
{
    structure(list(
        "statistic" = structure(x$F.0, names = x$termlabels),
        "permutations" = x$F.perm,
        "alternative" = "greater"),
              class = "permustats")
}

## no summary: test is a summary by itself

`permustats.protest` <-
    function(x, ...)
{
    structure(list(
        "statistic" = structure(x$t0, names = "r"),
        "permutations" = x$t,
        "alternative" = "greater"),
              class = "permustats")
}

## uses summary.procrustes

### the following do not return permutation data
`permustats.CCorA` <-
    function(x, ...)
{
    stop("no permutation data available")
}

`permustats.envfit` <-
    function(x, ...)
{
    stop("no permutation data available")
}

`permustats.factorfit` <-
    function(x, ...)
{
    stop("no permutation data available")
}

`permustats.vectorfit` <-
    function(x, ...)
{
    stop("no permutation data available")
}

`permustats.mso` <-
    function(x, ...)
{
    stop("no permutation data available")
}

`permustats.permutest.betadisper` <-
    function(x, ...)
{
    ntypes <- NCOL(x$perm)
    alt <- if (ntypes > 1) {
        c("greater", rep("two.sided", ntypes - 1))
    } else {
        "greater"
    }
    structure(list("statistic" = x$statistic,
                   "permutations" = x$perm,
                   "alternative" = alt),
              class = "permustats")
}

## no summary: print is sufficient

`permustats.anova.cca` <-
    function(x, ...)
{
    if (is.null(attr(x, "F.perm")))
        stop("no permutation data available")
    F.perm <- attr(x, "F.perm")
    k <- !is.na(x$F)
    F.0 <- x$F[k]
    structure(list(
       "statistic" = structure(F.0, names = rownames(x)[k]),
       "permutations" = F.perm,
       "alternative" = "greater"),
       class = "permustats")
}

## no summary: anova is a summary by itself
