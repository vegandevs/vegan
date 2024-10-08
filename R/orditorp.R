`orditorp` <-
    function (x, display, labels, choices = c(1, 2), priority, select,
              cex = 0.7, pcex, col = par("col"),
              pcol, pch = par("pch"), air = 1, ...)
{
    if (missing(pcex))
        pcex <- cex
    if (missing(pcol))
        pcol <- col
    if (inherits(x, "ordiplot"))
        display <- match.arg(display, names(x))
    sco <- scores(x, display = display, choices = choices, ...)
    kk <- complete.cases(sco)
    if (missing(labels))
        labels <- rownames(sco)
    if (missing(priority))
        priority <- rowSums((scale(sco)^2))
    if (!missing(select)) {
        sco <- .checkSelect(select, sco)
        labels <- .checkSelect(select, labels)
        priority <- .checkSelect(select, priority)
        kk <- .checkSelect(select, kk)
    }
    ## remove NA scores
    sco <- sco[kk,]
    priority <- priority[kk]
    labels <- labels[kk]
    w <- abs(strwidth(labels, cex = cex))/2 * air
    h <- abs(strheight(labels, cex = cex))/2 * air
    xx <- cbind(sco[, 1] - w, sco[, 1] + w, sco[, 2] - h, sco[, 2] +
                h)
    is.na(priority) <- w == 0
    ord <- rev(order(priority, na.last = FALSE))
    xx <- xx[ord, ]
    sco <- sco[ord, ]
    labels <- labels[ord]
    tt <- logical(nrow(xx))
    tt[1] <- TRUE
    for (i in 2:(nrow(xx) - sum(is.na(priority)))) {
        j <- 1:(i - 1)
        j <- j[tt[j]]
        tt[i] <- all(xx[i, 1] > xx[j, 2] | xx[j, 1] > xx[i, 2] |
                     xx[i, 3] > xx[j, 4] | xx[j, 3] > xx[i, 4])
    }
    if (sum(!tt)) {
        if (length(pch) > 1)
            pch <- (pch[ord])[!tt]
        if (length(pcex) > 1)
            pcex <- (pcex[ord])[!tt]
        if (length(pcol) > 1)
            pcol <- (pcol[ord])[!tt]
        ordiArgAbsorber(sco[!tt, , drop = FALSE], pch = pch, cex = pcex,
                        col = pcol, FUN = points, ...)
    }
    if (length(cex) > 1)
        cex <- (cex[ord])[tt]
    if (length(col) > 1)
        col <- (col[ord])[tt]
    ordiArgAbsorber(sco[tt, , drop = FALSE], labels = labels[tt], cex = cex,
                    col = col, FUN = text, ...)
    names(tt) <- labels
    if (inherits(x, "ordiplot")) {
        attr(x[[display]], "orditorp") <- tt[order(ord)]
        invisible(x)
    } else {
        invisible(tt[order(ord)])
    }
}
