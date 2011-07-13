"ordihull" <-
    function (ord, groups, display = "sites",
              draw = c("lines", "polygon", "none"),
              show.groups, label = FALSE, ...)
{
    draw <- match.arg(draw)
    pts <- scores(ord, display = display, ...)
    if (!missing(show.groups)) {
        take <- groups %in% show.groups
        pts <- pts[take, , drop = FALSE]
        groups <- groups[take]
    }
    out <- seq(along = groups)
    inds <- names(table(groups))
    res <- list()
    ## Remove NA scores
    kk <- complete.cases(pts)
    for (is in inds) {
        gr <- out[groups == is & kk]
        if (length(gr) > 1) {
            X <- pts[gr, ]
            hpts <- chull(X)
            hpts <- c(hpts, hpts[1])
            if (draw == "lines")
                ordiArgAbsorber(X[hpts, ], FUN = lines, ...)
            else if (draw == "polygon")
                ordiArgAbsorber(X[hpts,], FUN = polygon, ...)
            if (label && draw != "none") {
                cntr <- colMeans(X)
                ordiArgAbsorber(cntr[1], cntr[2], labels = is,
                                FUN = text, ...)
            }
            res[[is]] <- X[hpts,]
        }
    }
    class(res) <- "ordihull"
    invisible(res)
}
