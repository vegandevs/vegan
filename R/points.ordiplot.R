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
    ## draw arrows when requested, also for "species" etc. R warns
    ## loudly on (almost) zero-length arrows shorter than 1/1000
    ## inches: do not draw them.
    if (arrows) {
        seeit <- diff(par("usr")[1:2])/par("pin")[1]/1000
        zeroarr <- sqrt(rowSums(sco^2)) < seeit
        arrows(0, 0, sco[!zeroarr,1], sco[!zeroarr,2], length = length, ...)
    } else {
        points(sco, ...)
    }
    invisible(x)
}
