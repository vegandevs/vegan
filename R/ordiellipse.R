"ordiellipse" <-
    function (ord, groups, display = "sites", kind = c("sd", "se"),
              conf, draw = c("lines", "polygon"), w = weights(ord, display),
              show.groups, ...)
{
    ## Define Circle for an ellipse: taken from the 'car' package
    theta <- (0:51) * 2 * pi/51
    Circle <- cbind(cos(theta), sin(theta))
    weights.default <- function(object, ...) NULL
    kind <- match.arg(kind)
    draw <- match.arg(draw)
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
        if (length(gr) > 2) {
            X <- pts[gr, ]
            W <- w[gr]
            mat <- cov.wt(X, W)
            if (kind == "se")
                mat$cov <- mat$cov/mat$n.obs
            if (missing(conf))
                t <- 1
            else t <- sqrt(qchisq(conf, 2))
            xy <- t(mat$center + t * t(Circle %*% chol(mat$cov)))
            if (draw == "lines")
                ordiArgAbsorber(xy, FUN = lines, ...)
            else 
                ordiArgAbsorber(xy[, 1], xy[, 2], FUN = polygon, ...)
        }
    }
    invisible()
}
