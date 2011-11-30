## Ordimedian finds the spatial medians for groups. Spatial medians
## are L1 norms or statistics that minimize sum of distances of points
## from the statistic and 1d they are the medians. The current
## algorithm minimizes the L1 norm with optim and is pretty
## inefficient. Package ICSNP has a better algorithm (and we may steal
## it from them later).
`ordimedian` <-
    function(ord, groups, display = "sites", label = FALSE, ...)
{
    ## Sum of distances from the statistic
    medfun <-
        function(x, ord) sum(sqrt(rowSums(sweep(ord, 2, x)^2)),
                              na.rm = TRUE)
    ## derivative of medfun (if NULL, optim will use numerical
    ## differentiation)
    dmedfun <- function(x, ord) {
        up <- -sweep(ord, 2, x)
        dn <- sqrt(rowSums(sweep(ord, 2, x)^2))
        colSums(sweep(up, 1, dn, "/"))
    }
    #dmedfun <- NULL
    pts <- scores(ord, display = display, ...)
    inds <- names(table(groups))
    medians <- matrix(NA, nrow = length(inds), ncol = ncol(pts))
    rownames(medians) <- inds
    colnames(medians) <- colnames(pts)
    for (i in inds) {
        X <- pts[groups == i, , drop = FALSE]
        if (NROW(X) > 0)
            medians[i, ] <- optim(apply(X, 2, median, na.rm = TRUE),
                              fn = medfun, gr = dmedfun,
                              ord = X, method = "BFGS")$par
        if(label)
            ordiArgAbsorber(medians[i,1], medians[i,2], label = i,
                            FUN = text, ...)
    }
    invisible(medians)
}
