## Ordimedian finds the spatial medians for groups. Spatial medians
## are central locations that minimize sum of distances of points from
## the statistic and in 1D they are the medians. The current algorithm
## minimizes sum of distances with optim and is pretty
## inefficient. Package ICSNP has a better algorithm (and we may steal
## it from them later). Popular Weiszfeld's algorithm uses inverse
## distances as weights and iterates median as weighted average. Repeated
## application of the following will converge to spatial medians (and
## same results as this function):
##     w <- 1/sqrt(rowSums(X, 2, med)^2)
##     med <- colSums(w * X) / sum(w)
## but fails with infinite w when med is exactly on an data point.
`ordimedian` <-
    function(ord, groups, w, display = "sites", label = FALSE, ...)
{
    if (missing(w) || is.null(w))
        w <- rep_len(1, length(groups))
    ## Sum of distances from the statistic
    medfun <-
        function(x, ord, w) sum(w * sqrt(rowSums(sweep(ord, 2, x)^2)),
                              na.rm = TRUE)
    ## derivative of medfun (if NULL, optim will use numerical
    ## differentiation)
    dmedfun <- function(x, ord, w) {
        up <- -sweep(ord, 2, x)
        dn <- sqrt(rowSums(sweep(ord, 2, x)^2))
        colSums(w * sweep(up, 1, dn, "/"))
    }
    #dmedfun <- NULL
    pts <- scores(ord, display = display, ...)
    inds <- names(table(groups))
    medians <- matrix(NA, nrow = length(inds), ncol = ncol(pts))
    dimnames(medians) <- list(inds, colnames(pts))
    for (i in inds) {
        X <- pts[groups == i, , drop = FALSE]
        wgr <- w[groups == i]
        if (NROW(X) > 0)
            medians[i, ] <- optim(apply(X, 2, median, na.rm = TRUE),
                              fn = medfun, gr = dmedfun, w = wgr,
                              ord = X, method = "BFGS")$par
        if(label)
            ordiArgAbsorber(medians[i,1], medians[i,2], label = i,
                            FUN = text, ...)
    }
    invisible(medians)
}
