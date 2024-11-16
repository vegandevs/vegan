`points.ordiplot`  <-
    function (x, what, select, arrows = FALSE, length = 0.05,
              arr.mul, ...)
{
    sco <- scores(x, display = what)
    if (!missing(select))
        sco <- .checkSelect(select, sco)
    if (!missing(arr.mul)) {
        arrows <- TRUE
        sco <- sco * arr.mul
    } else {
        ## draw adjusted arrows automatically for biplot scores
        scoatt <- attr(sco, "score")
        if (!is.null(scoatt) && scoatt %in% c("biplot", "regression")) {
            arrows = TRUE
            sco <- sco * ordiArrowMul(sco)
        }
    }
    ## draw arrows when requested, also for "species" etc
    if (arrows) {
        arrows(0, 0, sco[,1], sco[,2], length = length, ...)
    } else {
        points(sco, ...)
    }
    invisible(x)
}
