plot.ordipointlabel <- function (x, ...)
{
    plot(x$points, pch = x$args$pch, cex = x$args$pcex, col = x$args$pcol,
        bg = x$args$pbg, asp = 1, ...)
    font <- attr(x$labels, "font")
    if (is.null(font))
        font <- par("font")
    text(x$labels, rownames(x$labels), cex = x$args$tcex, col = x$args$tcol,
         font = font, ...)
    psize <- par("din")
    if(any(abs(psize - x$dim)/x$dim > 0.1))
        message(gettextf(
            "original plot size was %.1f x %.1f, current is %.1f x %.1f",
            x$dim[1], x$dim[2], psize[1], psize[2]))
    invisible(x)
}
