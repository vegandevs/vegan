`points.procrustes` <-
    function(x, display = c("target","rotated"), truemean = FALSE, ...)
{
    display <- match.arg(display)
    X <- if (display == "target") x$X else x$Yrot
    if (truemean)
        X <- sweep(X, 2, x$xmean, "+")
    ordiArgAbsorber(X, FUN = points, ...)
    invisible()
}

`text.procrustes` <-
    function(x, display = c("target","rotated"), labels, truemean = FALSE, ...)
{
    display <- match.arg(display)
    X <- if (display == "target") x$X else x$Yrot
    if (truemean)
        X <- sweep(X, 2, x$xmean, "+")
    if (missing(labels))
        labels <- rownames(X)
    ordiArgAbsorber(X, labels = labels, FUN = text, ...)
    invisible()
}
