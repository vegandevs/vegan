panel.ordi <-
function(x, y, biplot, polygon, type = type, ...)
{
    tp <- trellis.par.get()
    sp <- tp$superpose.symbol
    ps <- tp$plot.symbol
    ## ordixyplot passes polygon of all points, but ordisplom does not
    if ("polygon" %in% type && !missing(polygon)) {
        ppar <- tp$plot.polygon
        lpolygon(polygon, col = ppar$col, border = ppar$border,
                 alpha = ppar$alpha, lty = ppar$lty, lwd = ppar$lwd, ...)
        inpol <- chull(x, y)
        par <- tp$superpose.polygon
        lpolygon(x[inpol], y[inpol], col = par$col[1L], border = par$border[1L],
                 alpha = par$alpha[1L], lty = par$lty[1L], lwd = par$lwd[1L], ...)
    }
    panel.xyplot(x, y, type = type,  ...)
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
        panel.superpose(x, y, panel.groups= "panel.ordiarrows", ...)
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
    col <- trellis.par.get("superpose.line")$col
    group <- list(...)$group.number
    if (is.null(group))
        col <- col[1]
    else
        col = rep(col, len=group)[group]
    panel.arrows(x[-n], y[-n], x[-1], y[-1], ends = ends, type = "open",
                 length = length, angle = angle, code = code,
                 col = col,
                 lwd = trellis.par.get("superpose.line")$lwd,
                 )
}
