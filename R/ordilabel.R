`ordilabel` <-
    function(x, display, labels, choices = c(1,2), priority, select,
             cex = 0.8, fill = "white", border = NULL, col = NULL,
             xpd = TRUE, ...)
{
    if (missing(display))
        display <- "sites"
    x <- scores(x, choices = choices, display = display, ...)
    if (missing(labels))
        labels <- rownames(x)
    if (!missing(select)) {
        x <- .checkSelect(select, x)
        labels <- .checkSelect(select, labels)
    }
    if (!missing(priority)) {
        if (!missing(select))
            priority <- priority[select]
        ord <- order(priority)
        x <- x[ord, ]
        labels <- labels[ord]
    }
    em <- strwidth("m", cex = cex, ...)
    ex <- strheight("x", cex = cex, ...)
    w <- (strwidth(labels, cex=cex,...) + em/1.5)/2
    h <- (strheight(labels, cex = cex, ...) + ex/1.5)/2
    if (is.null(col))
        if (!is.null(border))
            col <- border
        else
            col <- par("fg")
    for (i in 1:nrow(x)) {
        polygon(x[i,1] + c(-1,1,1,-1)*w[i], x[i,2] + c(-1,-1,1,1)*h[i],
                col = fill, border = border, xpd = xpd)
        text(x[i,1], x[i,2], labels = labels[i], cex = cex, col = col,
             xpd = xpd, ...)
    }
    invisible(x)
}

