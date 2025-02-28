`plot.ordisurf` <- function(x, what = c("contour","persp","gam"),
                            add = FALSE, bubble = FALSE, col = "red", cex = 1,
                            nlevels = 10, levels, labcex = 0.6,
                            lwd.cl = par("lwd"), ...) {
    what <- match.arg(what)
    y <- x$model$y
    x1 <- x$model$x1
    x2 <- x$model$x2
    X <- x$grid$x
    Y <- x$grid$y
    Z <- x$grid$z
    force(col)
    force(cex)
    if(what == "contour") {
        if(!add) {
            if(bubble) {
                if (is.numeric(bubble))
                    cex <- bubble
                cex <- (y -  min(y))/diff(range(y)) * (cex-0.4) + 0.4
            }
            plot(x1, x2, asp = 1, cex = cex, ...)
        }
        if (missing(levels))
            levels <- pretty(range(x$grid$z, finite = TRUE), nlevels)
        contour(X, Y, Z, col = col, add = TRUE,
                levels = levels, labcex = labcex,
                drawlabels = !is.null(labcex) && labcex > 0,
                lwd = lwd.cl)
    } else if(what == "persp") {
        persp(X, Y, Z, col = col, cex = cex, ...)
    } else {
        class(x) <- class(x)[-1]
        plot(x, ...) ##col = col, cex = cex, ...)
        class(x) <- c("ordisurf", class(x))
    }
    invisible(x)
}
