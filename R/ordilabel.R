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
    } else {
        ord <- seq_along(labels)
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
    col <- rep(col, length=nrow(x))[ord]
    if(!is.null(border))
        border <- rep(border, length=nrow(x))[ord]
    fill <- rep(fill, length=nrow(x))[ord]
    for (i in 1:nrow(x)) {
        ordiArgAbsorber(x[i,1] + c(-1,1,1,-1)*w[i], x[i,2] + c(-1,-1,1,1)*h[i],
                        col = fill[i], border = border[i], xpd = xpd,
                        FUN = polygon, ...)
        ordiArgAbsorber(x[i,1], x[i,2], labels = labels[i], cex = cex,
                        col = col[i], xpd = xpd, FUN = text, ...)
    }
    invisible(x)
}

