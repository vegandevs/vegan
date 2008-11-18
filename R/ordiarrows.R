"ordiarrows" <-
    function (ord, groups, levels, replicates, display = "sites",
              show.groups, startmark, ...)
{
    pts <- scores(ord, display = display, ...)
    npoints <- nrow(pts)
    if (missing(groups))
        groups <- gl(levels, replicates, npoints)
    if (!missing(show.groups)) {
        take <- groups %in% show.groups
        pts <- pts[take, , drop = FALSE]
        groups <- groups[take]
    }
    out <- seq(along = groups)
    inds <- names(table(groups))
    for (is in inds) {
        gr <- out[groups == is]
        if (length(gr) > 1) {
            X <- pts[gr, , drop = FALSE]
            X0 <- X[-nrow(X), , drop = FALSE]
            X1 <- X[-1, , drop = FALSE]
            nseg <- nrow(X0)
            if (!missing(startmark))
                points(X0[1,1], X0[1,2], pch=startmark, ...)
            if (nseg > 1)
                ordiArgAbsorber(X0[-nseg,1], X0[-nseg,2], X1[-nseg,1],
                                X1[-nseg,2], FUN = segments, ...)
            ordiArgAbsorber(X0[nseg, 1], X0[nseg, 2], X1[nseg, 1], X1[nseg, 2],
                            FUN = arrows, ...)
        }
    }
    invisible()
}
