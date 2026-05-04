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
    ## R warns on "zero-length" arrows shorter than 1/1000 inches in
    ## the current device: do not draw them, see help(arrows)
    if (arrows) { # optimize w.r.t. arrowheads
        seeit <- diff(par("usr")[1:2])/par("pin")[1]/1000
        zeroarr <- sqrt(rowSums(sco^2)) < seeit
        arrows(0, 0, sco[!zeroarr,1], sco[!zeroarr,2], length = length, ...)
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
