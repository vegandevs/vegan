`points.procrustes` <-
    function(x, display = c("target","rotated"), choices = c(1,2),
             truemean = FALSE, ...)
{
    display <- match.arg(display)
    X <- if (display == "target") x$X else x$Yrot
    X <- X[, choices, drop = FALSE]
    if (truemean)
        X <- sweep(X, 2, x$xmean[choices], "+")
    ordiArgAbsorber(X, FUN = points, ...)
    invisible()
}

`text.procrustes` <-
    function(x, display = c("target","rotated"), choices = c(1,2),
             labels, truemean = FALSE, ...)
{
    display <- match.arg(display)
    X <- if (display == "target") x$X else x$Yrot
    X <- X[, choices, drop = FALSE]
    if (truemean)
        X <- sweep(X, 2, x$xmean[choices], "+")
    if (missing(labels))
        labels <- rownames(X)
    ordiArgAbsorber(X, labels = labels, FUN = text, ...)
    invisible()
}
