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
    function(object, probs, ...)
{
    ## default cut levels for quantiles: these are two-sided
    if (missing(probs))
        probs <- switch(object$alternative,
                        "two.sided" = c(0.025, 0.5, 0.975),
                        "greater" = c(0.5, 0.95),
                        "less" = c(0.05, 0.5)) 
    sim <- t(object$permutations)
    object$means <- rowMeans(sim)
    sd <- apply(sim, 1, sd)
    object$z <-
        (object$statistic - object$means)/sd
    object$quantile <-
        apply(sim, 1, quantile, probs = probs, na.rm = TRUE)
    ## not (yet) P-values...
    class(object) <- "summary.permustats"
    object
}

`print.summary.permustats` <-
    function(x, ...)
{
    m <- cbind("statistic" = x$statistic,
               "z" = x$z,
               "mean" = x$means,
               t(x$quantile))
    printCoefmat(m, cs.ind = 3:ncol(m), ...)
    invisible(x)
}

### densityplot

`densityplot.permustats` <-
    function(x, data, xlab = "Permutations", ...)
{
    obs <- x$statistic
    sim <- rbind(x$statistic, as.matrix(x$permutations))
    nm <- names(obs)[col(sim)]
    densityplot( ~ as.vector(sim) | factor(nm, levels = unique(nm)),
                xlab = xlab,
                panel = function(x, ...) {
                    panel.densityplot(x, ...)
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

`qqmath.permustats` <-
    function(x, data, observed = TRUE, ylab = "Permutations", ...)
{
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

###
### specific methods to extract permustats
###

`permustats.anosim` <-
    function(x, ...)
{
    structure(list(
        "statistic" = structure(x$statistic, names="R"),
        "permutations" = x$perm,
        "alternative" = "greater"),
              class="permustats")
}

`permustats.adonis` <-
    function(x, ...)
{
    tab <- x$aov.tab
    k <- !is.na(tab$F.Model)
    structure(list(
        "statistic" = structure(tab$F.Model[k], names = rownames(tab)[k]),
        "permutations" = x$f.perms,
        "alternative" = "greater"),
              class="permustats")
}

`permustats.mantel` <-
    function(x, ...)
{
    structure(list(
        "statistic" = structure(x$statistic, names="r"),
        "permutations" = x$perm,
        "alternative" = "greater"),
              class="permustats")
}

`permustats.mrpp` <-
    function(x, ...)
{
    structure(list(
        "statistic" = structure(x$delta, names="delta"),
        "permutations" = x$boot.deltas,
        "alternative" = "less"),
              class="permustats")
}

`permustats.oecosimu` <-
    function(x, ...)
{
    structure(list(
        "statistic" = x$oecosimu$statistic,
        "permutations" = t(x$oecosimu$simulated),
        "alternative" = x$oecosimu$alternative),
              class="permustats")
}

`permustats.permutest.cca` <-
    function(x, ...)
{
    structure(list(
        "statistic" = structure(x$F.0, names = "F"),
        "permutations" = x$F.perm,
        "alternative" = "greater"),
              class="permustats")
}

`permustats.protest` <-
    function(x, ...)
{
    structure(list(
        "statistic" = structure(x$t0, names = "r"),
        "permutations" = x$t,
        "alternative" = "greater"),
              class="permustats")
}

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
    stop("no permutation data available")
}
