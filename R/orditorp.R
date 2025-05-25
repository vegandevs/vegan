`orditorp` <-
    function (x, display, labels, choices = c(1, 2), priority, select,
              cex = 0.7, pcex, col = par("col"),
              pcol, pch = par("pch"), air = 1, ...)
{
    if (missing(pcex))
        pcex <- cex
    if (missing(pcol))
        pcol <- col
    if (PIPE <- inherits(x, "ordiplot")) # PIPE returns x
        display <- match.arg(display, names(x))
    sco <- scores(x, display = display, choices = choices, ...)
    kk <- complete.cases(sco)
    if (!missing(select)) {
        names(kk) <- rownames(sco)
        sco <- .checkSelect(select, sco)
        kk <- .checkSelect(select, kk)
        if (!missing(priority) && length(priority) > NROW(sco))
            priority <- .checkSelect(select, priority)
    }
    if (missing(labels))
        labels <- rownames(sco)
    else
        rownames(sco) <- labels
    if (missing(priority))
        priority <- rowSums((scale(sco)^2))
    ## remove NA scores
    sco <- sco[kk,, drop = FALSE]
    ## Handle pathological cases of no data and one label
    if (identical(nrow(sco), 0L)) {
        if (PIPE)
            return(invisible(x))
        else
            stop("nothing to plot")
    }
    else if (identical(nrow(sco), 1L)) {
        ordiArgAbsorber(sco, labels = labels, cex = cex,
                        col = col, FUN = text, ...)
        tt <- structure(TRUE, names = labels)
        if (PIPE) {
            attr(x[[display]], "orditorp") <- tt
            return(invisible(x))
        } else
            return(invisible(tt))
    }
    priority <- priority[kk]
    labels <- labels[kk]
    w <- abs(strwidth(labels, cex = cex))/2 * air
    h <- abs(strheight(labels, cex = cex))/2 * air
    xx <- cbind(sco[, 1] - w, sco[, 1] + w, sco[, 2] - h, sco[, 2] +
                h)
    is.na(priority) <- w == 0
    ord <- rev(order(priority, na.last = FALSE))
    xx <- xx[ord,, drop = FALSE]
    sco <- sco[ord,, drop = FALSE]
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
    if (PIPE) {
        attr(x[[display]], "orditorp") <- tt[order(ord)]
        invisible(x)
    } else {
        invisible(tt[order(ord)])
    }
}
