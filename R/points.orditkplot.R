`points.orditkplot` <- function(x, pch = x$args$pch, cex = x$args$pcex,
                                col = x$args$pcol, bg = x$args$pbg, ...) {
    points(x$points, pch = pch, cex = cex, col = col, bg = bg, ...)
}

