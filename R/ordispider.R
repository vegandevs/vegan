"ordispider" <-
    function (ord, groups, display = "sites", w = weights(ord, display),
              show.groups, ...)
{
    weights.default <- function(object, ...) NULL
    if (inherits(ord, "cca") && missing(groups)) {
        lc <- scores(ord, display = "lc", ...)
        wa <- scores(ord, display = "wa", ...)
        ordiArgAbsorber(lc[, 1], lc[, 2], wa[, 1], wa[, 2],
                        FUN = segments, ...)
        return(invisible())
    }
    pts <- scores(ord, display = display, ...)
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
    out <- seq(along = groups)
    inds <- names(table(groups))
    for (is in inds) {
        gr <- out[groups == is]
        if (length(gr) > 1) {
            X <- pts[gr, ]
            W <- w[gr]
            ave <- apply(X, 2, weighted.mean, w = W)
            ordiArgAbsorber(ave[1], ave[2], X[, 1], X[, 2],
                            FUN = segments, ...)
        }
    }
    invisible()
}
