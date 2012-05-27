"points.metaMDS" <-
    function (x, display = c("sites", "species"),
              choices = c(1, 2), shrink = FALSE, select, ...)
{
    display <- match.arg(display)
    x <- scores(x, display = display, choices = choices, shrink = shrink)
    if (!missing(select))
        x <- .checkSelect(select, x)
    points(x, ...)
    invisible()
}
