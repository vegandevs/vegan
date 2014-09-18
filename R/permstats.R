### Functions to extract permutation statististic or null model
### results from various vegan objects.

## extract items as 'statistic' and 'simulations'. Specific methods
## towards the end of this file

`permstats` <-
    function(x, ...)
{
    UseMethod("permstats")
}

## something like str()
`print.permstats` <-
    function(x, ...)
{
    print(str(x))
    invisible(x)
}

### modelled after print.oecosimu (should perhaps have oecosimu() args
### like 'alternative'

`summary.permstats` <-
    function(object, ...)
{
    ## cut levels for quantiles: these are two-sided
    probs <- c(0.025, 0.5, 0.975)
    sim <- t(object$simulations)
    object$means <- rowMeans(sim)
    sd <- apply(sim, 1, sd)
    object$z <-
        (object$statistic - object$means)/sd
    object$quantile <-
        apply(sim, 1, quantile, probs = probs, na.rm = TRUE)
    ## not (yet) P-values...
    class(object) <- "summary.permstats"
    object
}

`print.summary.permstats` <-
    function(x, ...)
{
    m <- cbind("statistic" = x$statistic,
               "z" = x$z,
               "mean" = x$means,
               t(x$quantile))
    printCoefmat(m, cs.ind = 3:6, ...)
    invisible(x)
}

### specific methods

`permstats.anosim` <-
    function(x, ...)
{
    structure(list(
        "statistic" = x$statistic,
        "simulations" = x$perm),
              class="permstats")
}

`permstats.adonis` <-
    function(x, ...)
{
    stat <- x$aov.tab$F.Model
    stat <- stat[!is.na(stat)]
    structure(list(
        "statistic" = stat,
        "simulations" = x$f.perms),
              class="permstats")
}

`permstats.mantel` <-
    function(x, ...)
{
    structure(list(
        "statistic" = x$statistic,
        "simulations" = x$perm),
              class="permstats")
}

`permstats.mrpp` <-
    function(x, ...)
{
    structure(list(
        "statistic" = x$delta,
        "simulations" = x$boot.deltas),
              class="permstats")
}

`permstats.oecosimu` <-
    function(x, ...)
{
    structure(list(
        "statistic" = x$oecosimu$statistic,
        "simulations" = t(x$oecosimu$simulated)),
              class="permstats")
}

`permstats.permutest.cca` <-
    function(x, ...)
{
    structure(list(
        "statistic" = x$F.0,
        "simulations" = x$F.perm),
              class="permstats")
}

`permstats.protest` <-
    function(x, ...)
{
    structure(list(
        "statistic" = x$t.0,
        "simulations" = x$t),
              class="permstats")
}

### the following do not return permutation data
`permstats.CCorA` <-
    function(x, ...)
{
    stop("no permutation data available")
}

`permstats.envfit` <-
    function(x, ...)
{
    stop("no permutation data available")
}

`permstats.factorfit` <-
    function(x, ...)
{
    stop("no permutation data available")
}

`permstats.vectorfit` <-
    function(x, ...)
{
    stop("no permutation data available")
}

`permstats.mso` <-
    function(x, ...)
{
    stop("no permutation data available")
}

`permstats.permutest.betadisper` <-
    function(x, ...)
{
    stop("no permutation data available")
}
