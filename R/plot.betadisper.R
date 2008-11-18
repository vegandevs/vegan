`plot.betadisper` <- function(x, axes = c(1,2), cex = 0.7, hull = TRUE,
                              ylab, xlab, main, sub, ...)
{
    localAxis <- function(..., col, bg, pch, cex, lty, lwd) axis(...)
    localBox <- function(..., col, bg, pch, cex, lty, lwd) box(...)
    localTitle <- function(..., col, bg, pch, cex, lty, lwd) title(...)
    if(missing(main))
        main <- deparse(substitute(x))
    if(missing(sub))
        sub <- paste("method = \"", attr(x, "method"), "\"", sep = "")
    if(missing(xlab))
        xlab <- paste("PCoA", axes[1])
    if(missing(ylab))
        ylab <- paste("PCoA", axes[2])
    g <- scores(x, choices = axes)
    plot(g$sites, asp = 1, type = "n", axes = FALSE, ann = FALSE, ...)
    for(i in levels(x$group)) {
        j <- which(levels(x$group) == i)
        segments(g$centroids[j, axes[1]],
                 g$centroids[j, axes[2]],
                 g$sites[x$group == i, axes[1]],
                 g$sites[x$group == i, axes[2]], col = "blue", ...)
        if(hull) {
            ch <- chull(g$sites[x$group == i, axes])
            ch <- c(ch, ch[1])
            lines(x$vectors[x$group == i, axes][ch, ],
                       col = "black", lty = "dashed", ...)
        }
    }
    points(g$centroids, pch = 16, cex = 1, col = "red", ...)
    points(g$sites, pch = as.numeric(x$group),
           cex = cex, ...)
    localTitle(main = main, xlab = xlab, ylab = ylab, sub = sub, ...)
    localAxis(1, ...)
    localAxis(2, ...)
    localBox(...)
    class(g) <- "ordiplot"
    invisible(g)
}
