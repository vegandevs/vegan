`points.orditkplot` <- function(x, pch = x$args$pch, cex = x$args$pcex,
                                col = x$args$pcol, bg = x$args$pbg, ...)
{
    .Deprecated(msg = "function was moved to CRAN package vegan3d")
    points(x$points, pch = pch, cex = cex, col = col, bg = bg, ...)
}

