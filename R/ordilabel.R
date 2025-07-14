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
    cex <- rep(cex, length.out = nrow(sco))
    if (!missing(priority)) {
        ord <- order(priority)
        sco <- sco[ord, ]
        labels <- labels[ord]
        cex <- cex[ord]
    } else {
        ord <- seq_along(labels)
    }
    ## allow variable cex
    w <- numeric(nrow(sco))
    h <- numeric(nrow(sco))
    for(i in seq_along(cex)) {
        em <- strwidth("m", cex = cex[i], ...)
        ex <- strheight("x", cex = cex[i], ...)
        w[i] <- (strwidth(labels[i], cex=cex[i],...) + em/1.5)/2
        h[i] <- (strheight(labels[i], cex = cex[i], ...) + ex/1.5)/2
    }
    if (is.null(col))
        col <- par("fg")
    col <- rep(col, length=nrow(sco))[ord]
    if (is.null(border))
        border <- col
    border <- rep(border, length=nrow(sco))[ord]
    fill <- rep(fill, length=nrow(sco))[ord]
    dev.hold()
    on.exit(dev.flush())
    for (i in 1:nrow(sco)) {
        ordiArgAbsorber(sco[i,1] + c(-1,1,1,-1)*w[i],
                        sco[i,2] + c(-1,-1,1,1)*h[i],
                        col = fill[i], border = border[i], xpd = xpd,
                        FUN = polygon, ...)
        ordiArgAbsorber(sco[i,1], sco[i,2], labels = labels[i], cex = cex[i],
                        col = col[i], xpd = xpd, FUN = text, ...)
    }
    invisible(x)
}

