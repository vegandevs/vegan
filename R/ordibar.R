### draws crossed error bars for classes in ordination. These are
### oblique to axis because so are th clouds of the points and their
### standard errors and confidence regions. The bars are principal
### axes of corresponding ellipse (as drawn in ordiellipse), and found
### as principal components of the associate covariance matrix. The
### function is modelled after ordiellipse.
`ordibar` <-
    function (ord, groups, display = "sites", kind = c("sd", "se"),
              conf,  w = weights(ord, display), col = 1,
              show.groups, label = FALSE, lwd = NULL, length = 0,  ...)
{
    weights.default <- function(object, ...) NULL
    kind <- match.arg(kind)
    draw <- TRUE
    pts <- scores(ord, display = display, ...)
    ## ordibar only works with 2D data (2 columns)
    pts <- as.matrix(pts)
    if (ncol(pts) > 2)
        pts <- pts[ , 1:2, drop = FALSE]
    if (ncol(pts) < 2)
        stop("needs two dimensions")
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
    if (label) {
        cntrs <- matrix(NA, nrow=length(inds), ncol=2)
        rownames(cntrs) <- inds
    }
    col <- rep(col, length = length(inds))
    names(col) <- inds
    res <- list()
    ## Remove NA scores
    kk <- complete.cases(pts) & !is.na(groups)
    for (is in inds) {
        gr <- out[groups == is & kk]
        if (length(gr)) {
            X <- pts[gr, , drop = FALSE]
            W <- w[gr]
            mat <- cov.wt(X, W)
            if (mat$n.obs == 1)
                mat$cov[] <- 0
            if (kind == "se")
                mat$cov <- mat$cov * sum(mat$wt^2)
            if (missing(conf))
                t <- 1
            else t <- sqrt(qchisq(conf, 2))
            if (mat$n.obs > 1) {
                eig <- eigen(mat$cov, symmetric = TRUE)
                v <- sweep(eig$vectors, 2, sqrt(eig$values), "*") * t
                cnt <- mat$center
                ordiArgAbsorber(v[1,] + cnt[1], v[2,] + cnt[2],
                                -v[1,] + cnt[1], -v[2,] + cnt[2],
                                col = col[is], lwd = lwd,
                                length = length/2, angle = 90, code = 3,
                                FUN = arrows, ...)
            }
            if (label) {
                cntrs[is,] <- mat$center
            }
            mat$scale <- t
            res[[is]] <- mat
        }
    }
    if (label) {
        ordiArgAbsorber(cntrs, col = par("fg"), border = col,
                        FUN = ordilabel, ...)
    }
    class(res) <- "ordibar"
    invisible(res)
}
