`text.ordiplot`  <-
    function (x, what, labels, select, arrows = FALSE, length = 0.05, ...)
{
    sco <- scores(x, display = what)
    if (!missing(labels))
        rownames(sco) <- labels
    if (!missing(select))
        sco <- .checkSelect(select, sco)
    scoatt <- attr(sco, "score")
    if (!is.null(scoatt) && scoatt %in% c("biplot", "regression")) {
        arrows = TRUE
        sco <- sco * ordiArrowMul(sco)
    }
    if (arrows) {
        arrows(0, 0, sco[,1], sco[,2], length = length, ...)
        sco <- ordiArrowTextXY(sco, rownames(sco), rescale = FALSE, ...)
    }
    text(sco, labels = rownames(sco), ...)
    invisible(x)
}
