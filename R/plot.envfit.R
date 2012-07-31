`plot.envfit` <-
    function (x, choices = c(1, 2), labels, arrow.mul, at = c(0, 0),
              axis = FALSE, p.max = NULL, col = "blue", bg, add = TRUE, ...)
{
    formals(arrows) <- c(formals(arrows), alist(... = ))
    ## get labels
    labs <- list("v" = rownames(x$vectors$arrows),
                 "f" = rownames(x$factors$centroids))
    ## Change labels if user so wishes
    if (!missing(labels)) {
        ## input list of "vectors" and/or "factors"
        if (is.list(labels)) {
            if (!is.null(labs$v) && !is.null(labels$vectors))
                labs$v <- labels$vectors
            if (!is.null(labs$f) && !is.null(labels$factors))
                labs$f <- labels$factors
        } else {
            ## input vector: either vectors or factors must be NULL,
            ## and the existing set of labels is replaced
            if (!is.null(labs$v) && !is.null(labs$f))
                stop("needs a list with both 'vectors' and 'factors' labels")
            if (!is.null(labs$v))
                labs$v <- labels
            else
                labs$f <- labels
        }
    }
    vect <- NULL
    if (!is.null(p.max)) {
        if (!is.null(x$vectors)) {
            take <- x$vectors$pvals <= p.max
            x$vectors$arrows <- x$vectors$arrows[take, , drop = FALSE]
            labs$v <- labs$v[take]
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
            labs$f <- labs$f[take]
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
        plot.new() ## needed for string widths and heights
        if(!is.null(vect)) {
            ## compute axis limits allowing space for labels
            sw <- strwidth(labs$v, ...)
            sh <- strheight(labs$v, ...)
            xlim <- range(at[1], vtext[,1] + sw, vtext[,1] - sw)
            ylim <- range(at[2], vtext[,2] + sh, vtext[,2] - sh)
            if(!is.null(x$factors)) {
                ## if factors, also need to consider them
                sw <- strwidth(labs$f, ...)
                sh <- strheight(labs$f, ...)
                xlim <- range(xlim, x$factors$centroids[, choices[1]] + sw,
                              x$factors$centroids[, choices[1]] - sw)
                ylim <- range(ylim, x$factors$centroids[, choices[2]] + sh,
                              x$factors$centroids[, choices[2]] - sh)
            }
            ## these plotting calls will prob. generate warnings
            ## because of passing ... everywhere. localFoo needed?
            plot.window(xlim = xlim, ylim = ylim, asp = 1, ...)
            axis(side = 1, ...)
            axis(side = 2, ...)
            box(...)
            alabs <- colnames(vect)
            title(..., ylab = alabs[2], xlab = alabs[1])
        } else if (!is.null(x$factors)) {
            sw <- strwidth(labs$f, ...)
            sh <- strheight(labs$f, ...)
            xlim <- range(at[1], x$factors$centroids[, choices[1]] + sw,
                          x$factors$centroids[, choices[1]] - sw)
            ylim <- range(at[2], x$factors$centroids[, choices[2]] + sh,
                          x$factors$centroids[, choices[2]] - sh)
            ## these plotting calls will prob. generate warnings
            ## because of passing ... everywhere. localFoo needed?
            plot.window(xlim = xlim, ylim = ylim, asp = 1, ...)
            axis(side = 1, ...)
            axis(side = 2, ...)
            box(...)
            alabs <- colnames(x$factors$centroids[, choices, drop = FALSE])
            title(..., ylab = alabs[2], xlab = alabs[1])
        } else stop("Nothing to plot")
    }
    if (!is.null(vect)) {
        arrows(at[1], at[2], vect[, 1], vect[, 2], len = 0.05,
               col = col)
        if (missing(bg))
            text(vtext, labs$v, col = col, ...)
        else
            ordilabel(vtext, labels = labs$v, col = col, fill = bg, ...)
    }
    if (!is.null(x$factors)) {
        if (missing(bg))
            text(x$factors$centroids[, choices, drop = FALSE],
                 labs$f, col = col, ...)
        else
            ordilabel(x$factors$centroids[, choices, drop = FALSE],
                      labels = labs$f, col = col, fill = bg, ...)
    }
    if (axis && !is.null(vect)) {
        axis(3, at = ax + at[1], labels = c(maxarr, 0, maxarr),
             col = col)
        axis(4, at = ax + at[2], labels = c(maxarr, 0, maxarr),
             col = col)
    }
    invisible()
}
