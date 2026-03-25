`text.ordiplot`  <-
    function (x, what, labels, select, optimize = FALSE, arrows = FALSE,
              length = 0.05, arr.mul, bg, ...)
{
    sco <- scores(x, display = what)
    if (!missing(select))
        sco <- .checkSelect(select, sco)
    if (!missing(labels))
        rownames(sco) <- labels
    if (!missing(arr.mul)) {
        arrows <- TRUE
        sco <- sco * arr.mul
    } else {
        scoatt <- attr(sco, "score")
        if (!is.null(scoatt) && scoatt %in% c("biplot", "regression")) {
            arrows <- TRUE
            sco <- sco * ordiArrowMul(sco)
        }
    }
    if (arrows) { # optimize w.r.t. arrowheads
        arrows(0, 0, sco[,1], sco[,2], length = length, ...)
        if (!optimize)
            sco <- ordiArrowTextXY(sco, rownames(sco), rescale = FALSE, ...)
    }
    if (optimize) { # draw no points at arrowheads
        if (missing(bg))
            ordipointlabel(sco, display = what, add = TRUE, points = !arrows,
                           ...)
        else
            ordipointlabel(sco, display = what, bg = bg, add = TRUE,
                           points = !arrows, ...)
    } else if (missing(bg))
        text(sco, labels = rownames(sco), ...)
    else
        ordilabel(sco, labels = rownames(sco), fill = bg, ...)
    invisible(x)
}
