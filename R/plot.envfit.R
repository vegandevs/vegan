`plot.envfit` <-
    function (x, choices = c(1, 2), arrow.mul, at = c(0, 0), 
              axis = FALSE, p.max = NULL, col = "blue", add = TRUE, ...) 
{
    formals(arrows) <- c(formals(arrows), alist(... = ))
    vect <- NULL
    if (!is.null(p.max)) {
        if (!is.null(x$vectors)) {
            take <- x$vectors$pvals <= p.max
            x$vectors$arrows <- x$vectors$arrows[take, , drop = FALSE]
            x$vectors$r <- x$vectors$r[take]
            if (nrow(x$vectors$arrows) == 0) 
                x$vectors <- vect <- NULL
        }
        if (!is.null(x$factors)) {
            tmp <- x$factors$pvals <= p.max
            nam <- names(tmp)[tmp]
            take <- x$factors$var.id %in% nam
            x$factors$centroids <- x$factors$centroids[take, 
                                                       , drop = FALSE]
            if (nrow(x$factors$centroids) == 0) 
                x$factors <- NULL
        }
    }
    if (!is.null(x$vectors)) {
        vect <- sqrt(x$vectors$r) * x$vectors$arrows[, choices, drop = FALSE]
        if (missing(arrow.mul)) {
            if(!add)
                arrow.mul <- 1
            else 
                arrow.mul <- ordiArrowMul(vect, at = at)
        }
        if (axis) {
            maxarr <- round(sqrt(max(x$vectors$r)), 1)
            ax <- -c(-1, 0, 1) * arrow.mul * maxarr
        }
        vect <- arrow.mul * vect
        vtext <- sweep(1.1 * vect, 2, at, "+")
        vect <- sweep(vect, 2, at, "+")
    }
    if (!add) {
        xlim <- range(at[1], vect[,1], x$factors$centroids[,1])
        ylim <- range(at[2], vect[,2], x$factors$centroids[,2])
        if (!is.null(vect)) 
            plot(vect, xlim = xlim, ylim = ylim, asp = 1, type = "n", 
                 ...)
        else if (!is.null(x$factors)) 
            plot(x$factors$centroids[, choices, drop = FALSE], 
                 asp = 1, xlim = xlim, ylim = ylim, type = "n", 
                 ...)
        else stop("Nothing to plot")
    }
    if (!is.null(vect)) {
        arrows(at[1], at[2], vect[, 1], vect[, 2], len = 0.05, 
               col = col)
        text(vtext, rownames(x$vectors$arrows), col = col, ...)
    }
    if (!is.null(x$factors)) {
        text(x$factors$centroids[, choices, drop = FALSE], rownames(x$factors$centroids), 
             col = col, ...)
    }
    if (axis && !is.null(vect)) {
        axis(3, at = ax + at[1], labels = c(maxarr, 0, maxarr), 
             col = col)
        axis(4, at = ax + at[2], labels = c(maxarr, 0, maxarr), 
             col = col)
    }
    invisible()
}
