`plot.ordisurf` <- function(x, what = c("contour","surface", "persp","gam"),
                            add = FALSE, bubble = FALSE, col = "red",
                            alpha = 1, cex = 1, nlevels = 10, levels,
                            labcex = 0.6, lwd.cl = par("lwd"), ...)
{
    what <- match.arg(what)
    y <- x$model$y
    x1 <- x$model$x1
    x2 <- x$model$x2
    X <- x$grid$x
    Y <- x$grid$y
    Z <- x$grid$z
    force(col)
    force(cex)
    if(what %in% c("contour", "surface")) {
        if (what == "surface" && length(col) < 2) # default image colors
            col <- hcl.colors(12, "YlOrRd", rev = TRUE, alpha = alpha)
        if (!add)
            plot(X, Y, asp = 1, type="n", ...)
        if (what == "surface")
            image(X, Y, Z, add = TRUE, col = col)
        if (missing(levels))
            levels <- pretty(range(x$grid$z, finite = TRUE), nlevels)
        contour(X, Y, Z, col = if (what == "surface") par("fg") else col[1],
                add = TRUE, levels = levels, labcex = labcex,
                drawlabels = !is.null(labcex) && labcex > 0,
                lwd = lwd.cl)
        if(bubble && !add) {
            if (is.numeric(bubble))
                cex <- bubble
            cex <- (y -  min(y))/diff(range(y)) * (cex-0.4) + 0.4
            points(x1, x2, cex = cex, ...)
        }
    } else if(what == "persp") {
        persp(X, Y, Z, col = col, cex = cex, ...)
    } else {
        class(x) <- class(x)[-1]
        plot(x, asp =1, ...) ##col = col, cex = cex, ...)
        class(x) <- c("ordisurf", class(x))
    }
    invisible(x)
}
