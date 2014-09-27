### Functions to extract permutation statististic or null model
### results from various vegan objects.

## extract items as 'statistic' and 'simulations'. Specific methods
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
    function(object, ...)
{
    ## cut levels for quantiles: these are two-sided
    probs <- switch(object$alternative,
                    "two.sided" = c(0.025, 0.5, 0.975),
                    "greater" = c(0.5, 0.95),
                    "less" = c(0.05, 0.5)) 
    sim <- t(object$simulations)
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
    function(x, xlab = "Permutations", ...)
{
    obs <- x$statistic
    sim <- rbind(x$statistic, as.matrix(x$simulations))
    nm <- names(obs)[col(sim)]
    densityplot( ~ as.vector(sim) | factor(nm, levels = unique(nm)),
                xlab = xlab,
                panel = function(x, ...) {
                    panel.densityplot(x, ...)
                    panel.abline(v = obs[panel.number()], ...)
                },
                ...)
}


### specific methods

`permustats.anosim` <-
    function(x, ...)
{
    structure(list(
        "statistic" = structure(x$statistic, names="R"),
        "simulations" = x$perm,
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
        "simulations" = x$f.perms,
        "alternative" = "greater"),
              class="permustats")
}

`permustats.mantel` <-
    function(x, ...)
{
    structure(list(
        "statistic" = structure(x$statistic, names="r"),
        "simulations" = x$perm,
        "alternative" = "greater"),
              class="permustats")
}

`permustats.mrpp` <-
    function(x, ...)
{
    structure(list(
        "statistic" = structure(x$delta, names="delta"),
        "simulations" = x$boot.deltas,
        "alternative" = "less"),
              class="permustats")
}

`permustats.oecosimu` <-
    function(x, ...)
{
    structure(list(
        "statistic" = x$oecosimu$statistic,
        "simulations" = t(x$oecosimu$simulated),
        "alternative" = x$oecosimu$alternative),
              class="permustats")
}

`permustats.permutest.cca` <-
    function(x, ...)
{
    structure(list(
        "statistic" = structure(x$F.0, names = "F"),
        "simulations" = x$F.perm,
        "alternative" = "greater"),
              class="permustats")
}

`permustats.protest` <-
    function(x, ...)
{
    structure(list(
        "statistic" = structure(x$t0, names = "r"),
        "simulations" = x$t,
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
