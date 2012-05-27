"points.ordiplot" <-
    function (x, what, select, ...)
{
    x <- scores(x, what)
    if (!missing(select))
        x <- .checkSelect(select, x)
    points(x, ...)
    invisible()
}
