"points.ordiplot" <-
    function (x, what, select, ...) 
{
    x <- scores(x, what)
    if (!missing(select))
        x <- x[select,,drop=FALSE]
    points(x, ...)
    invisible()
}
