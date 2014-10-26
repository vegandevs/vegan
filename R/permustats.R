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

`summary.permustats` <- function(object, interval = 0.95, ...) {
    nalt <- length(object$alternative)
    nstat <- length(object$statistic)
    ## Replicate alternative to length of statistic
    if ((nalt < nstat) && identical(nalt, 1L)) {
        object$alternative <- rep(object$alternative, length.out = nstat)
    }
    TAB <- c("two.sided", "greater", "less")
    compint <- (1 - interval) / 2
    PROBS <- list(two.sided = c(compint, 0.5, interval + compint),
                  greater = c(NA, 0.5, interval),
                  less = c(1 - interval, 0.5, NA))
    alt <- match(object$alternative, TAB)
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
    ## not (yet) P-values...
    class(object) <- "summary.permustats"
    object
}

`print.summary.permustats` <- function(x, ...) {
    m <- cbind("statistic" = x$statistic,
               "z" = x$z,
               "mean" = x$means,
               x$quantile)
    cat("\n")
    printCoefmat(m, tst.ind = 1:ncol(m), na.print = "", ...)
    writeLines(strwrap(paste0("(Interval (Upper - Lower) = ", x$interval, ")", sep = ""),
                       initial = "\n"))
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

`permustats.ordiareatest` <-
    function(x, ...)
{
    structure(list(
        "statistic" = x$areas,
        "permutations" = t(x$permutations),
        "alternative" = "less"),
              class = "permustats")
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
    ntypes <- NCOL(x$perm)
    alt <- if (ntypes > 1) {
        c("greater", rep("two.sided", ntypes - 1))
    } else {
        "greater"
    }
    structure(list("statistic" = x$statistic,
                   "permutations" = x$perm,
                   "alternative" = alt),
              class ="permustats")
}
