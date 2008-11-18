"points.decorana" <-
    function (x, display=c("sites", "species"), choices=1:2, origin = TRUE,
              select, ...)
{
    display <- match.arg(display)
    x <- scores(x, display = display, choices = choices, origin = origin, ...)
    if (!missing(select))
        x <- x[select,,drop=FALSE]
    ordiArgAbsorber(x, FUN = points, ...)
    invisible()
}

