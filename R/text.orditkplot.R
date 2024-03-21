`text.orditkplot` <- function(x, cex = x$args$tcex, col = x$args$tcol,
                              font = attr(x$labels, "font"), ...)
{
    .Deprecated(msg = "function was moved to CRAN package vegan3d")
    if (is.null(font)) {
        font <- par("font")
    }
    text(x$labels, labels = rownames(x$labels), cex = cex, col = col,
         font = font, ...)
}

