`ordihull` <-
    function (ord, groups, display = "sites",
              draw = c("lines", "polygon", "none"),
              col = NULL, show.groups, label = FALSE, ...)
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
    if (label)
        cntrs <- names <- NULL
    ## Remove NA scores
    kk <- complete.cases(pts)
    for (is in inds) {
        gr <- out[groups == is & kk]
        if (length(gr) > 1) {
            X <- pts[gr, ]
            hpts <- chull(X)
            hpts <- c(hpts, hpts[1])
            if (draw == "lines")
                ordiArgAbsorber(X[hpts, ], FUN = lines,
                                col = if(is.null(col)) par("fg") else col, ...)
            else if (draw == "polygon")
                ordiArgAbsorber(X[hpts,], FUN = polygon, col = col, ...)
            if (label && draw != "none") {
                cntrs <- rbind(cntrs, colMeans(X[hpts[-1],, drop = FALSE]))
                names <- c(names, is)
            }
            res[[is]] <- X[hpts,]
        }
    }
    if (label && draw != "none") {
        if (draw == "lines")
            ordiArgAbsorber(cntrs[,1], cntrs[,2], labels = names,
                            col = col, FUN = text, ...)
        else
            ordiArgAbsorber(cntrs, labels = names, col = NULL,
                            FUN = ordilabel, ...)
    }
    class(res) <- "ordihull"
    invisible(res)
}
