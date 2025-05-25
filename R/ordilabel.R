`ordilabel` <-
    function(x, display, labels, choices = c(1,2), priority, select,
             cex = 0.8, fill = "white", border = NULL, col = NULL,
             xpd = TRUE, ...)
{
    if (missing(display))
        display <- "sites"
    sco <- scores(x, choices = choices, display = display, ...)
    if (!missing(select)) {
        sco <- .checkSelect(select, sco)
        if (!missing(priority) && length(priority) > NROW(sco))
            priority <- .checkSelect(select, priority)
    }
    if (missing(labels))
        labels <- rownames(sco)
    if (!missing(priority)) {
        ord <- order(priority)
        sco <- sco[ord, ]
        labels <- labels[ord]
    } else {
        ord <- seq_along(labels)
    }
    em <- strwidth("m", cex = cex, ...)
    ex <- strheight("x", cex = cex, ...)
    w <- (strwidth(labels, cex=cex,...) + em/1.5)/2
    h <- (strheight(labels, cex = cex, ...) + ex/1.5)/2
    if (is.null(col))
        col <- par("fg")
    col <- rep(col, length=nrow(sco))[ord]
    if (is.null(border))
        border <- col
    border <- rep(border, length=nrow(sco))[ord]
    fill <- rep(fill, length=nrow(sco))[ord]
    dev.hold()
    for (i in 1:nrow(sco)) {
        ordiArgAbsorber(sco[i,1] + c(-1,1,1,-1)*w[i],
                        sco[i,2] + c(-1,-1,1,1)*h[i],
                        col = fill[i], border = border[i], xpd = xpd,
                        FUN = polygon, ...)
        ordiArgAbsorber(sco[i,1], sco[i,2], labels = labels[i], cex = cex,
                        col = col[i], xpd = xpd, FUN = text, ...)
    }
    dev.flush()
    invisible(x)
}

