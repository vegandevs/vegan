`plot.specaccum` <-
    function(x, add = FALSE, random = FALSE, ci = 2,
             ci.type = c("bar","line","polygon"), col = par("fg"), lty = 1,
             ci.col = col, ci.lty = 1, ci.length = 0,  xlab, ylab = x$method,
             ylim, xvar = c("sites", "individuals", "effort"), ...)
{
    if(random && !(x$method %in% c("random", "collector")))
        stop("random = TRUE can be used only with method='random'")
    if(x$method == "collector")
        random <- TRUE
    xvar <- match.arg(xvar)
    ## adjust weights to number of sites
    if (random && !is.null(x$weights) && xvar == "sites") {
        n <- length(x$effort)
        adj <- n/x$effort[n]
    } else {
        adj <- 1
    }
    xaxvar <- x[[xvar]]
    if (missing(xlab))
        xlab <- paste(toupper(substring(xvar, 1, 1)),
                              substring(xvar, 2), sep="")
    if (random)
        ci <- FALSE
    ci.type <- match.arg(ci.type)
    if (!add) {
        if (missing(ylim))
            if (random)
                ylim <- c(min(1, x$perm), max(x$perm, na.rm = TRUE))
            else
                ylim <- c(min(1, x$richness),
                          max(x$richness, x$richness + ci*x$sd, na.rm = TRUE))
        plot(xaxvar, x$richness, xlab=xlab, ylab=ylab, ylim=ylim,
             type="n", ...)
    }
    if (!is.null(x$sd) && ci)
        switch(ci.type,
               bar = {k <- x$sd > sqrt(.Machine$double.eps)
                      arrows(xaxvar[k], (x$richness - ci*x$sd)[k], xaxvar[k],
                             (x$richness + ci*x$sd)[k], col=ci.col, lty=ci.lty,
                             length = ci.length, angle = 90, code = 3, ...)},
               line = matlines(xaxvar, x$richness + t(rbind(-ci,ci) %*% x$sd),
                 col=ci.col, lty=ci.lty, ...),
               polygon = polygon(c(xaxvar, rev(xaxvar)),
                 c(x$richness - ci*x$sd, rev(x$richness + ci*x$sd)), col=ci.col,
                 lty=ci.lty,  ...)
               )
    if (random) {
        if (is.null(x$weights)) {
            for(i in seq_len(NCOL(x$perm)))
                lines(xaxvar, x$perm[,i], col=col, ...)
        } else {
            for(i in seq_len(NCOL(x$perm)))
                lines(x$weights[,i]*adj, x$perm[,i], col=col, ...)
        }
    } else
        lines(xaxvar, x$richness, col=col, lty = lty, ...)
    invisible()
}

`lines.specaccum` <-
    function(x, ...)
{
    plot(x, add = TRUE, ...)
}
