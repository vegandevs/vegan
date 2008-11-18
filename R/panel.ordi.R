`panel.ordi` <-
    function(x, y, biplot, type = type,  ...)
{
    panel.xyplot(x, y, type = type,  ...)
    tp <- trellis.par.get()
    if ("biplot" %in% type && !is.null(biplot$arrows)) {
        panel.arrows(0, 0, biplot$arrows[,1], biplot$arrows[,2],
                     col=tp$superpose.symbol$col, ...)
    }
    if ("biplot" %in% type && !is.null(biplot$centres)) {
        panel.xyplot(biplot$centres[,1], biplot$centres[,2],
                     col = tp$plot.symbol$col, 
                     pch = "+", cex = 3*tp$plot.symbol$cex, lwd=2,
                     ...)
    }
    panel.abline(h=0, lty = 3)
    panel.abline(v=0, lty = 3)
}
