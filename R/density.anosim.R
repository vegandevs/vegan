### density & densityplot methods for vegan functions returning
### statistics from permuted/simulated data. These are modelled after
### density.oecosimu and densityplot.oecosimu (which are in their
### separate files).

## anosim

`density.anosim` <-
    function(x, ...)
{
    out <- density(x$perm, ...)
    out$call <- match.call()
    out$call[[1]] <- as.name("density")
    out
}

## adonis can return a matrix of terms, hence we also have densityplot()

`density.adonis` <-
    function(x, ...)
{
    cols <- ncol(x$f.perms)
    if (cols > 1)
        warning("'density' is meaningful only with one term, you have ", cols)
    out <- density(x$f.perms, ...)
    out$call <- match.call()
    out$call[[1]] <- as.name("density")
    out
}

`densityplot.adonis` <-
    function(x, data, xlab = "Null", ...)
{
    require(lattice) || stop("requires package 'lattice'")
    sim <- x$f.perms
    obs <- x$aov.tab$F.Model
    nm <- rownames(x$aov.tab)[col(sim)]
    densityplot( ~ as.vector(sim) | factor(nm, levels = unique(nm)),
                xlab = xlab,
                panel = function(x, ...) {
                    panel.densityplot(x, ...)
                    panel.abline(v = obs[panel.number()], ...)
                },
                ...)
}

## mantel

`density.mantel` <-
    function(x, ...)
{
    out <- density(x$perm, ...)
    out$call <- match.call()
    out$call[[1]] <- as.name("density")
    out
}

## mrpp

`density.mrpp` <-
    function(x, ...)
{
    out <- density(x$boot.deltas, ...)
    out$call <- match.call()
    out$call[[1]] <- as.name("density")
    out
}

## anova.cca does not return permutation results, but permutest.cca
## does. However, permutest.cca always finds only one statistic. Full
## tables anova.cca are found by repeated calls to permutest.cca.

`density.permutest.cca` <-
    function(x, ...)
{
    out <- density(x$F.perm, ...)
    out$call <- match.call()
    out$call[[1]] <- as.name("density")
    out
}

## protest

`density.protest` <-
    function(x, ...)
{
    out <- density(x$t, ...)
    out$call <- match.call()
    out$call[[1]] <- as.name("density")
    out
}
