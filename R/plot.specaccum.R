"plot.specaccum" <-
    function(x, add = FALSE, ci = 2, ci.type = c("bar","line","polygon"), 
             col = par("fg"), ci.col = col, ci.lty = 1, xlab = "Sites",
             ylab = x$method, ylim, ...)
{
    ci.type <- match.arg(ci.type)
    if (!add) {
        if (missing(ylim))
            ylim<- c(1, max(x$richness, x$richness + ci*x$sd))
        plot(x$sites, x$richness, xlab=xlab, ylab=ylab, ylim=ylim,
             type="n", ...)
    }
    if (!is.null(x$sd) && ci)
        switch(ci.type,
               bar = segments(x$sites, x$richness - ci*x$sd, x$sites,
                  x$richness + ci*x$sd, col=ci.col, lty=ci.lty, ...),
               line = matlines(x$sites, x$richness + t(rbind(-ci,ci) %*% x$sd),
                 col=ci.col, lty=ci.lty, ...),
               polygon = polygon(c(x$sites, rev(x$sites)),
                 c(x$richness - ci*x$sd, rev(x$richness + ci*x$sd)), col=ci.col,
                 lty=ci.lty,  ...)
               )
    lines(x$sites, x$richness,col=col, ...)
    invisible()
}
