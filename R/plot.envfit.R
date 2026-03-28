`plot.envfit` <-
    function (x, choices = c(1, 2), labels, arrow.mul, at = c(0, 0),
              axis = FALSE, p.max = NULL, r2.min = NULL, col = "blue", bg,
              optimize = FALSE, cex = 1, add = TRUE, ...)
{
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
            if (!is.null(x$vectors) && !is.null(x$factors))
                stop("needs a list with both 'vectors' and 'factors' labels")
            ## need to handle the case where both sets of labels are NULL
            ## such as when used with the default interface and single x
            if (!is.null(x$factors))
                labs$f <- labels
            else
                labs$v <- labels
        }
    }
    vect <- NULL
    if (!is.null(p.max) || !is.null(r2.min)) {
        if (!is.null(x$vectors)) {
            take <- x$vectors$pvals <= min(p.max, 1) &
                x$vectors$r >= max(r2.min, 0)
            x$vectors$arrows <- x$vectors$arrows[take, , drop = FALSE]
            labs$v <- labs$v[take]
            x$vectors$r <- x$vectors$r[take]
            if (nrow(x$vectors$arrows) == 0)
                x$vectors <- vect <- NULL
        }
        if (!is.null(x$factors)) {
            tmp <- x$factors$pvals <= min(p.max, 1) &
                x$factors$r >= max(r2.min, 0)
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
    }
    ## vectors may need scaling (arrow.mul) and shifting (at)
    if (!is.null(vect)) {
        if (missing(arrow.mul)) {
            if (add) {
                arrow.mul <- ordiArrowMul(vect, at = at)
            } else {
                if (is.null(x$factors))
                    arrow.mul <- 1
                else {
                    plot.window(xlim = range(x$factors$centroids[,1], at[1]),
                                ylim = range(x$factors$centroids[,2], at[2]),
                                asp = 1, ...)
                    arrow.mul <- ordiArrowMul(vect, at = at)
                }
            }
        }
        vect <- arrow.mul * vect
        if (!missing(at))
            vect <- sweep(vect, 2, at, "+")
    }
    ## Create a new plot instead of adding to an existing one
    if (!add) {
        ordiplot(rbind(vect, x$factors$centroids), display = "sites",
                 type = "n", ...)
    }
    ## add elements to the existing plot
    ## arrows for vectors
    if (!is.null(vect))
        arrows(at[1], at[2], vect[, 1], vect[, 2], length = 0.05,
               col = col)
    ## text labels
    pch <- "x"
    if (!optimize && !is.null(vect)) {
        vect <- ordiArrowTextXY(vect, labs$v, rescale = FALSE, at = at, cex=cex,
                                ...)
    } else {
        pch <- c(rep("x", NROW(x$factors$centroids)),
                 rep("", NROW(vect)))
    }
    text.ordiplot(rbind(x$factors$centroids, vect), "sites", col = col,
                  bg = bg, optimize = optimize, pch = pch, cex = cex, ...)
    if (axis && !is.null(vect)) {
        maxr <- max(sqrt(abs(x$vectors$r)))
        ax <- pretty(c(-maxr, maxr))
        axis(3, at = arrow.mul * ax + at[1], labels = ax)
        axis(4, at = arrow.mul * ax + at[2], labels = ax)
    }
    invisible()
}
