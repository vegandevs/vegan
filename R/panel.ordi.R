panel.ordi <-
function(x, y, biplot, type = type,  ...)
{
    panel.xyplot(x, y, type = type,  ...)
    tp <- trellis.par.get()
    sp <- tp$superpose.symbol
    ps <- tp$plot.symbol
    if ("biplot" %in% type && !is.null(biplot$arrows)) {
        panel.arrows(0, 0, biplot$arrows[,2], biplot$arrows[,1],
                     col=sp$col, ...)
    }
    if ("biplot" %in% type && !is.null(biplot$centres)) {
        panel.xyplot(biplot$centres[,2], biplot$centres[,1],
                     col = ps$col,
                     pch = "+", cex = 3*ps$cex, lwd=2,
                     ...)
    }
    if ("arrows" %in% type) {
        panel.superpose(x, y, panel.groups= "panel.ordiarrows",  ...)
    }
    panel.abline(h=0, lty = 3)
    panel.abline(v=0, lty = 3)
}

## needed for "arrows" %in% type
panel.ordiarrows <-
function(x, y, subscripts,
         ends = "last", type = "open", length = 0.25, angle = 30, code = 2,
         ...)
{
    n <- length(x)
    panel.arrows(x[-n], y[-n], x[-1], y[-1], ends = ends, type = "open",
                 length = length, angle = angle, code = code,
                 col = trellis.par.get("superpose.line")$col,
                 lwd = trellis.par.get("superpose.line")$lwd,
                 )
}
