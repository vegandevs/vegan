"ordiellipse" <-
    function (ord, groups, display = "sites", kind = c("sd", "se"),
              conf, draw = c("lines", "polygon", "none"),
              w = weights(ord, display),
              show.groups, label = FALSE,  ...)
{
    weights.default <- function(object, ...) NULL
    kind <- match.arg(kind)
    draw <- match.arg(draw)
    pts <- scores(ord, display = display, ...)
    ## ordiellipse only works with 2D data (2 columns)
    pts <- as.matrix(pts)
    if (ncol(pts) > 2)
        pts <- pts[ , 1:2, drop = FALSE]
    if (ncol(pts) < 2)
        stop("ordiellipse needs two dimensions")
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
    res <- list()
    if (label)
        cntrs <- names <- NULL
    for (is in inds) {
        gr <- out[groups == is]
        if (length(gr) > 2) {
            X <- pts[gr, ]
            W <- w[gr]
            mat <- cov.wt(X, W)
            if (kind == "se")
                mat$cov <- mat$cov/mat$n.obs
            if (missing(conf))
                t <- 1
            else t <- sqrt(qchisq(conf, 2))
            xy <- veganCovEllipse(mat$cov, mat$center, t)
            if (draw == "lines")
                ordiArgAbsorber(xy, FUN = lines, ...)
            else if (draw == "polygon") 
                ordiArgAbsorber(xy[, 1], xy[, 2], FUN = polygon, ...)
            if (label && draw != "none") {
                cntrs <- rbind(cntrs, mat$center)
                names <- c(names, is)
            }
            mat$scale <- t
            res[[is]] <- mat
        }
    }
    if (label && draw != "none") {
        if (draw == "lines")
            ordiArgAbsorber(cntrs[,1], cntrs[,2], labels=names, 
                            FUN = text, ...)
        else
            ordiArgAbsorber(cntrs, labels = names, FUN = ordilabel, ...)
    }
    class(res) <- "ordiellipse"
    invisible(res)
}
