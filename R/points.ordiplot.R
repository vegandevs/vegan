`points.ordiplot`  <-
    function (x, what, select, ...)
{
    sco <- scores(x, what)
    if (!missing(select))
        sco <- .checkSelect(select, sco)
    points(sco, ...)
    invisible(x)
}
