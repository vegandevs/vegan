`biplot.CCorA` <-
    function(x, xlabs, which = 1:2, ...)
{
    plottable <- x$Mat.ranks[which] > 1
    if (!all(plottable)) {
        warning("plot ", paste(which(!plottable), collapse=","),
                " not drawn because it only has one dimension")
        which <- which[plottable]
    }
    if (prod(par("mfrow")) < length(which)) {
        op <- par(mfrow=c(1,2))
        on.exit(par(op))
    }
    if (missing(xlabs))
        xlabs <- rownames(x$Cy)
    else if (!is.null(xlabs) && is.na(xlabs))
        xlabs <- rep(NA, nrow(x$Cy))
    else if (is.null(xlabs))
        xlabs <- 1:nrow(x$Cy)
    if (any(which == 1))
        biplot(x$Cy, x$AA, xlabs = xlabs, ...)
    if (any(which == 2))
        biplot(x$Cx, x$BB, xlabs = xlabs, ...)
    invisible()
}
