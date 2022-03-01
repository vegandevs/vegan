`points.ordiplot`  <-
    function (x, what, select, arrows = FALSE, ...)
{
    sco <- scores(x, what)
    if (!missing(select))
        sco <- .checkSelect(select, sco)
    ## draw adjusted arrows automatically for biplot scores
    scoatt <- attr(sco, "score")
    if (!is.null(scoatt) && scoatt %in% c("biplot", "regression")) {
        arrows = TRUE
        sco <- sco * ordiArrowMul(sco)
    }
    ## draw arrows when requested, also for "species" etc
    if (arrows) {
        arrows(0, 0, sco[,1], sco[,2], ...)
    } else {
        points(sco, ...)
    }
    invisible(x)
}
