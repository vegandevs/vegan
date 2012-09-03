### density & densityplot methods for vegan functions returning
### statistics from permuted/simulated data. These are modelled after
### density.oecosimu and densityplot.oecosimu (which are in their
### separate files).

## anosim

`density.anosim` <-
    function(x, ...)
{
    obs <- x$statistic
    ## Put observed statistic among permutations
    out <- density(c(obs, x$perm), ...)
    out$call <- match.call()
    out$observed <- obs
    out$call[[1]] <- as.name("density")
    class(out) <- c("vegandensity", class(out))
    out
}

## adonis can return a matrix of terms, hence we also have densityplot()

`density.adonis` <-
    function(x, ...)
{
    cols <- ncol(x$f.perms)
    if (cols > 1)
        warning("'density' is meaningful only with one term, you have ", cols)
    obs <- x$aov.tab$F.Model
    obs <- obs[!is.na(obs)]
    out <- density(c(obs, x$f.perms), ...)
    out$observed <- obs
    out$call <- match.call()
    out$call[[1]] <- as.name("density")
    class(out) <- c("vegandensity", class(out))
    out
}

`densityplot.adonis` <-
    function(x, data, xlab = "Null", ...)
{
    require(lattice) || stop("requires package 'lattice'")
    obs <- x$aov.tab$F.Model
    obs <- obs[!is.na(obs)]
    sim <- rbind(obs, x$f.perms)
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
    obs <- x$statistic
    out <- density(c(obs, x$perm), ...)
    out$observed <- obs
    out$call <- match.call()
    out$call[[1]] <- as.name("density")
    class(out) <- c("vegandensity", class(out))
    out
}

## mrpp

`density.mrpp` <-
    function(x, ...)
{
    obs <- x$delta
    out <- density(c(obs, x$boot.deltas), ...)
    out$observed <- obs
    out$call <- match.call()
    out$call[[1]] <- as.name("density")
    class(out) <- c("vegandensity", class(out))
    out
}

## anova.cca does not return permutation results, but permutest.cca
## does. However, permutest.cca always finds only one statistic. Full
## tables anova.cca are found by repeated calls to permutest.cca.

`density.permutest.cca` <-
    function(x, ...)
{
    obs <- x$F.0
    out <- density(c(obs, x$F.perm), ...)
    out$observed <- obs
    out$call <- match.call()
    out$call[[1]] <- as.name("density")
    class(out) <- c("vegandensity", class(out))
    out
}

## protest

`density.protest` <-
    function(x, ...)
{
    obs <- x$t0
    out <- density(c(obs, x$t), ...)
    out$observed <- obs
    out$call <- match.call()
    out$call[[1]] <- as.name("density")
    class(out) <- c("vegandensity", class(out))
    out
}

#### plot method: the following copies stats::plot.density() code but
#### adds one new argument to draw abline(v=...) for the observed
#### statistic

`plot.vegandensity` <-
    function (x, main = NULL, xlab = NULL, ylab = "Density", type = "l", 
    zero.line = TRUE, obs.line = TRUE, ...) 
{
    if (is.null(xlab)) 
        xlab <- paste("N =", x$n, "  Bandwidth =", formatC(x$bw))
    if (is.null(main)) 
        main <- deparse(x$call)
    ## change obs.line to col=2 (red) if it was logical TRUE
    if (isTRUE(obs.line))
        obs.line <- 2
    plot.default(x, main = main, xlab = xlab, ylab = ylab, type = type,
                 ...)
    if (zero.line) 
        abline(h = 0, lwd = 0.1, col = "gray")
    if (is.character(obs.line) || obs.line)
        abline(v = x$observed, col = obs.line)
    invisible(NULL)
}
