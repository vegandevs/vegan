`points.ordiplot`  <-
    function (x, what, select, arrows = FALSE, ...)
{
    sco <- scores(x, what)
    if (!missing(select))
        sco <- .checkSelect(select, sco)
    if (arrows) {
        sco <- sco * ordiArrowMul(sco)
        arrows(0, 0, sco[,1], sco[,2], ...)
    } else {
        points(sco, ...)
    }
    invisible(x)
}
