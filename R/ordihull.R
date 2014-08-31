`ordihull` <-
    function (ord, groups, display = "sites",
              draw = c("lines", "polygon", "none"),
              col = NULL, alpha = 127, show.groups, label = FALSE, ...)
{
    draw <- match.arg(draw)
    ## Internal function to find the polygon centre
    polycentre <- function(x) {
        n <- nrow(x)
        if (n < 4) 
            return(colMeans(x[-n, ]))
        xy <- x[-n, 1] * x[-1, 2] - x[-1, 1] * x[-n, 2]
        A <- sum(xy)/2
        xc <- sum((x[-n, 1] + x[-1, 1]) * xy)/A/6
        yc <- sum((x[-n, 2] + x[-1, 2]) * xy)/A/6
        c(xc, yc)
    }
    ## Make semitransparent fill colour
    if (draw == "polygon" && !is.null(col))
        col <- rgb(t(col2rgb(col)), alpha = alpha, maxColorValue = 255)
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
    kk <- complete.cases(pts) & !is.na(groups)
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
                cntrs <- rbind(cntrs, polycentre(X[hpts,]))
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
