`ordisegments` <-
    function (ord, groups, levels, replicates, order.by, display = "sites",
              col = 1, show.groups, label = FALSE, ...)
{
    pts <- scores(ord, display = display, ...)
    npoints <- nrow(pts)
    if (missing(groups))
        groups <- gl(levels, replicates, npoints)
    if (!missing(order.by)) {
        if (length(order.by) != nrow(pts))
            stop(gettextf("the length of order.by (%d) does not match the number of points (%d)",
                 length(order.by), nrow(pts)))
        ord <- order(order.by)
        pts <- pts[ord,]
        groups <- groups[ord]
    }
    if (!missing(show.groups)) {
        take <- groups %in% show.groups
        pts <- pts[take, , drop = FALSE]
        groups <- groups[take]
    }
    out <- seq(along = groups)
    inds <- names(table(groups))
    if (is.factor(col))
        col <- as.numeric(col)
    col <- rep(col, length=length(inds))
    names(col) <- inds
    ends <- names <- NULL
    for (is in inds) {
        gr <- out[groups == is]
        if (length(gr) > 1) {
            X <- pts[gr, , drop = FALSE]
            X0 <- X[-nrow(X), , drop = FALSE]
            X1 <- X[-1, , drop = FALSE]
            ordiArgAbsorber(X0[, 1], X0[, 2], X1[, 1], X1[, 2],
                            col = col[is], FUN = segments, ...)
            if (label) {
                ends <- rbind(ends, X[c(1, nrow(X)), ])
                names <- c(names, is, is)
            }
        }
    }
    if (label)
        ordiArgAbsorber(ends, labels = names, border = col, col = par("fg"),
                        FUN = ordilabel, ...)
    invisible()
}
