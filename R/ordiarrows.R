`ordiarrows` <-
    function (ord, groups, levels, replicates, order.by, 
              display = "sites", col = 1, show.groups, startmark,
              label = FALSE, length = 0.1, ...)
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
    starts <- names <- NULL
    for (is in inds) {
        gr <- out[groups == is]
        if (length(gr) > 1) {
            X <- pts[gr, , drop = FALSE]
            X0 <- X[-nrow(X), , drop = FALSE]
            X1 <- X[-1, , drop = FALSE]
            nseg <- nrow(X0)
            if (!missing(startmark))
                points(X0[1,1], X0[1,2], pch=startmark, col = col[is], ...)
            if (label) {
                starts <- rbind(starts, X0[1,])
                names <- c(names, is)
            }
            if (nseg > 1)
                ordiArgAbsorber(X0[-nseg,1], X0[-nseg,2], X1[-nseg,1],
                                X1[-nseg,2], col = col[is],
                                FUN = segments, ...)
            ordiArgAbsorber(X0[nseg, 1], X0[nseg, 2], X1[nseg, 1],
                            X1[nseg, 2], col = col[is], length = length,
                            FUN = arrows, ...)
        }
    }
    if (label)
        ordiArgAbsorber(starts, labels = names, border = col, col = par("fg"),
                        FUN = ordilabel, ...)
    invisible()
}
