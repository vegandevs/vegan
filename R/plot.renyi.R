### Simple base::plot that draws lines of all sites in one plot

`plot.renyi` <-
    function(x, ...)
{
    ## one site: make into 1-row matrix
    if (!is.data.frame(x))
        x <- t(x)
    divs <- colnames(x)
    ndivs <- length(divs)
    matplot(seq_len(ndivs), t(x), axes=FALSE, type="l",
            col=palette(), lty=rep(1:6, each=length(palette())),
            xlab = "alpha", ylab = "Diversity", ...)
    axis(1, at = seq_len(ndivs), labels=divs)
    axis(2)
    box()
    if(nrow(x) > 1)
        legend("topright", rownames(x), lty=rep(1:6, each=length(palette())),
                                                col = palette())
}

### Simple base::plot that draws accumulation by sites for all scales
### in one plot

`plot.renyiaccum` <-
    function(x, ...)
{
    ## Take only mean
    if ("permutation" %in% names(dimnames(x)))
        x[,,1] <- apply(x, 2, rowMeans)
    x <- x[,,1]
    matplot(seq_len(nrow(x)), x, type="l", lty=1, xlab = "Sites",
            ylab = "Diversity", ...)
    legend("topleft", colnames(x), lty=1, col=palette())
}

### Add lines for given statistic ("mean", "stdev", "min"...,
### "Collector") for every scale

`lines.renyiaccum` <-
    function(x, what, ...)
{
    what <- match.arg(what, dimnames(x)[[3]])
    matlines(seq_len(nrow(x)), x[,, what], col = palette(), ...)
}

### Simple base::plot for poolaccum() and estaccumR(). Draws only the
### mean value of permutations. For empirical CI, use
### ggvegan::autoplot.

`plot.poolaccum` <-         # also works with estaccumR
    function(x, ...)
{
    x <- x$means
    matplot(x[,1], x[,-1], type="l", lty=1, xlab = "Sites", ylab = "Richness",
            ...)
    legend("bottomright", colnames(x[,-1]), col = palette(), lty=1, ...)
}
