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
    probs <- c(0.025, 0.5, 0.975)
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
    printCoefmat(m, cs.ind = 3:6, ...)
    invisible(x)
}

### specific methods

`permustats.anosim` <-
    function(x, ...)
{
    structure(list(
        "statistic" = x$statistic,
        "simulations" = x$perm),
              class="permustats")
}

`permustats.adonis` <-
    function(x, ...)
{
    stat <- x$aov.tab$F.Model
    stat <- stat[!is.na(stat)]
    structure(list(
        "statistic" = stat,
        "simulations" = x$f.perms),
              class="permustats")
}

`permustats.mantel` <-
    function(x, ...)
{
    structure(list(
        "statistic" = x$statistic,
        "simulations" = x$perm),
              class="permustats")
}

`permustats.mrpp` <-
    function(x, ...)
{
    structure(list(
        "statistic" = x$delta,
        "simulations" = x$boot.deltas),
              class="permustats")
}

`permustats.oecosimu` <-
    function(x, ...)
{
    structure(list(
        "statistic" = x$oecosimu$statistic,
        "simulations" = t(x$oecosimu$simulated)),
              class="permustats")
}

`permustats.permutest.cca` <-
    function(x, ...)
{
    structure(list(
        "statistic" = x$F.0,
        "simulations" = x$F.perm),
              class="permustats")
}

`permustats.protest` <-
    function(x, ...)
{
    structure(list(
        "statistic" = x$t.0,
        "simulations" = x$t),
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
