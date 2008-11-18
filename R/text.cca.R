`text.cca` <-
    function (x, display = "sites", labels, choices = c(1, 2), scaling = 2, 
              arrow.mul, head.arrow = 0.05, select, ...) 
{
    formals(arrows) <- c(formals(arrows), alist(... = ))
    if (length(display) > 1) 
        stop("Only one `display' item can be added in one command.")
    pts <- scores(x, choices = choices, display = display, scaling = scaling)
    if (!missing(labels))
        rownames(pts) <- labels
    if (!missing(select)) 
        pts <- pts[select, , drop = FALSE]
    if (display == "cn") {
        cnam <- rownames(pts)
        text(pts, labels = cnam, ...)
        pts <- scores(x, choices = choices, display = "bp", scaling = scaling)
        bnam <- rownames(pts)
        pts <- pts[!(bnam %in% cnam), , drop = FALSE]
        if (nrow(pts) == 0) 
            return(invisible())
        else display <- "bp"
    }
    if (display == "bp") {
        if (missing(arrow.mul)) {
            arrow.mul <- ordiArrowMul(pts)
        }
        pts <- pts * arrow.mul
        arrows(0, 0, pts[, 1], pts[, 2], length = head.arrow, 
               ...)
        pts <- pts * 1.1
        axis(3, at = c(-arrow.mul, 0, arrow.mul), labels = rep("", 
                                                  3))
        axis(4, at = c(-arrow.mul, 0, arrow.mul), labels = c(-1, 
                                                  0, 1))
    }
    text(pts, labels = rownames(pts), ...)
    invisible()
}
