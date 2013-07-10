plot.ordipointlabel <- function (x, ...)
{
    plot(x$points, pch = x$args$pch, cex = x$args$pcex, col = x$args$pcol,
        bg = x$args$pbg, asp = 1, ...)
    font <- attr(x$labels, "font")
    if (is.null(font))
        font <- par("font")
    text(x$labels, rownames(x$labels), cex = x$args$tcex, col = x$args$tcol,
         font = font, ...)
    invisible(x)
}
