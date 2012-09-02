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
    out$observed <- x$statistic
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
    out <- density(x$f.perms, ...)
    out$observed <- x$aov.tab$F.Model
    out$call <- match.call()
    out$call[[1]] <- as.name("density")
    class(out) <- c("vegandensity", class(out))
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
    out$observed <- x$statistic
    out$call <- match.call()
    out$call[[1]] <- as.name("density")
    class(out) <- c("vegandensity", class(out))
    out
}

## mrpp

`density.mrpp` <-
    function(x, ...)
{
    out <- density(x$boot.deltas, ...)
    out$observed <- x$delta
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
    out <- density(x$F.perm, ...)
    out$observed <- x$F.0
    out$call <- match.call()
    out$call[[1]] <- as.name("density")
    class(out) <- c("vegandensity", class(out))
    out
}

## protest

`density.protest` <-
    function(x, ...)
{
    out <- density(x$t, ...)
    out$observed <- x$t0
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
    zero.line = TRUE, obs.line = FALSE, ...) 
{
    if (is.null(xlab)) 
        xlab <- paste("N =", x$n, "  Bandwidth =", formatC(x$bw))
    if (is.null(main)) 
        main <- deparse(x$call)
    ## adjust xlim of obs.line if needed
    if (is.character(obs.line) || obs.line) {
        xlim <- range(c(x$x, x$observed), na.rm = TRUE)
        ## change obs.line to col=2 (red) if it was logical TRUE
        if (isTRUE(obs.line))
            obs.line <- 2
    } else {
        xlim <- NULL
    }
    ## check for explicit xlim in the call and use it if specified
    if(!is.null(match.call(expand.dots = FALSE)$...$xlim))
        plot.default(x, main = main, xlab = xlab, ylab = ylab, type = type,
                     ...)
    else
        plot.default(x, main = main, xlab = xlab, ylab = ylab, type = type,
                     xlim = xlim, ...)
    if (zero.line) 
        abline(h = 0, lwd = 0.1, col = "gray")
    if (is.character(obs.line) || obs.line)
        abline(v = x$observed, col = obs.line)
    invisible(NULL)
}
