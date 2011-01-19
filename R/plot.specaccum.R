`plot.specaccum` <-
    function(x, add = FALSE, ci = 2, ci.type = c("bar","line","polygon"), 
             col = par("fg"), ci.col = col, ci.lty = 1, xlab,
             ylab = x$method, ylim, xvar = c("sites", "individuals"), ...)
{
    xvar <- match.arg(xvar)
    xaxvar <- x[[xvar]]
    if (missing(xlab))
        xlab <- paste(toupper(substring(xvar, 1, 1)),
                              substring(xvar, 2), sep="")
    ci.type <- match.arg(ci.type)
    if (!add) {
        if (missing(ylim))
            ylim <- c(1, max(x$richness, x$richness + ci*x$sd))
        plot(xaxvar, x$richness, xlab=xlab, ylab=ylab, ylim=ylim,
             type="n", ...)
    }
    if (!is.null(x$sd) && ci)
        switch(ci.type,
               bar = segments(xaxvar, x$richness - ci*x$sd, xaxvar,
                  x$richness + ci*x$sd, col=ci.col, lty=ci.lty, ...),
               line = matlines(xaxvar, x$richness + t(rbind(-ci,ci) %*% x$sd),
                 col=ci.col, lty=ci.lty, ...),
               polygon = polygon(c(xaxvar, rev(xaxvar)),
                 c(x$richness - ci*x$sd, rev(x$richness + ci*x$sd)), col=ci.col,
                 lty=ci.lty,  ...)
               )
    lines(xaxvar, x$richness,col=col, ...)
    invisible()
}
