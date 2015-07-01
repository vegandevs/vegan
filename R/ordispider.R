`ordispider` <-
    function (ord, groups, display = "sites", w = weights(ord, display),
              spiders = c("centroid", "median"),
              show.groups, label = FALSE, ...)
{
    weights.default <- function(object, ...) NULL
    spiders <- match.arg(spiders)
    if (inherits(ord, "cca") && missing(groups)) {
        lc <- scores(ord, display = "lc", ...)
        wa <- scores(ord, display = "wa", ...)
        ordiArgAbsorber(lc[, 1], lc[, 2], wa[, 1], wa[, 2],
                        FUN = segments, ...)
        class(lc) <- "ordispider"
        return(invisible(lc))
    }
    pts <- scores(ord, display = display, ...)
    ## spids stores pointwise centroids to be returned invisibly
    ## (transposed here so that filling is easier, but back-transposed
    ## when returned).
    spids <- t(array(NA, dim=dim(pts), dimnames = dimnames(pts)))
    ## ordihull: draw lines from centre to the points in the hull
    if (inherits(ord, "ordihull"))
        groups <- attr(pts, "hulls")
    w <- eval(w)
    if (length(w) == 1)
        w <- rep(1, nrow(pts))
    if (is.null(w))
        w <- rep(1, nrow(pts))
    if (!missing(show.groups)) {
        take <- groups %in% show.groups
        pts <- pts[take, , drop = FALSE]
        groups <- groups[take]
        w <- w[take]
    }
    if (spiders == "median" && sd(w) > sqrt(.Machine$double.eps))
        warning("weights are ignored with 'median' spiders")
    out <- seq(along = groups)
    inds <- names(table(groups))
    if (label) 
    cntrs <- names <- NULL
    ## 'kk' removes NA scores and NA groups
    kk <- complete.cases(pts) & !is.na(groups)
    for (is in inds) {
        gr <- out[groups == is & kk]
        if (length(gr)) {
            X <- pts[gr, , drop = FALSE]
            W <- w[gr]
            if (length(gr) > 1) {
                ave <- switch(spiders,
                              "centroid" = apply(X, 2, weighted.mean, w = W),
                              "median" = ordimedian(X, rep(1, nrow(X))))
                ordiArgAbsorber(ave[1], ave[2], X[, 1], X[, 2],
                                FUN = segments, ...)
            } else {
                ave <- X
            }
            spids[,gr] <- ave
            if (label) {
                cntrs <- rbind(cntrs, ave)
                names <- c(names, is)
            }
        }
    }
    if (label) 
        ordiArgAbsorber(cntrs, label = names, FUN = ordilabel, ...)
    spids <- t(spids)
    class(spids) <- "ordispider"
    invisible(spids)
}
