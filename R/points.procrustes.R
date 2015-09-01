`points.procrustes` <-
    function(x, display = c("target","rotated"), ...)
{
    display <- match.arg(display)
    x <- if (display == "target") x$X else x$Yrot
    points(x, ...)
    invisible()
}

`text.procrustes` <-
    function(x, display = c("target","rotated"), labels, ...)
{
    display <- match.arg(display)
    x <- if (display == "target") x$X else x$Yrot
    if (missing(labels))
        labels <- rownames(x)
    text(x, labels = labels, ...)
    invisible()
}
