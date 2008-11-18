"points.procrustes" <-
    function(x, display = c("target","rotated"), ...)
{
    display <- match.arg(display)
    x <- if (display == "target") x$X else x$Yrot
    points(x, ...)
    invisible()
}
