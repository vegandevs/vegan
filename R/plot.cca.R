`plot.cca` <- function (x, choices = c(1, 2), display = c("sp", "wa", "cn"),
                        scaling = "species", type, xlim, ylim, const,
                        correlation = FALSE, hill = FALSE, cex = 0.7, ...) {
    TYPES <- c("text", "points", "none")
    ## take care that bp arrows are also returned if only cn given
    if (any(display %in% c("c","cn")))
        display <- c(display, "bp")
    g <- scores(x, choices, display, scaling, const, correlation = correlation,
                hill = hill, tidy = FALSE, droplist = FALSE)
    if (length(g) == 0 || all(is.na(g)))
      stop("nothing to plot: requested scores do not exist")
    if (!is.list(g)) # never! but doesn't harm either
        g <- list(default = g)
    ## Take care that there are names
    for (i in seq_along(g)) {
        if (length(dim(g[[i]])) > 1)
            rownames(g[[i]]) <- rownames(g[[i]], do.NULL = FALSE,
                                         prefix = substr(names(g)[i], 1, 3))
    }
    if (!is.null(g$centroids)) {
        if (is.null(g$biplot)) # should never be null if g$centroids exist
            g$biplot <- scores(x, choices, "bp", scaling)
        bipnam <- rownames(g$biplot)
        cntnam <- rownames(g$centroids)
        g$biplot <- g$biplot[!(bipnam %in% cntnam), , drop = FALSE]
        if (nrow(g$biplot) == 0)
            g$biplot <- NULL
    }
    if (missing(type)) {
        nitlimit <- 80
        nit <- max(nrow(g$spe), nrow(g$sit), nrow(g$con), nrow(g$cen), 0)
        if (nit > nitlimit)
            type <- "points"
        else type <- "text"
    }
    else type <- match.arg(type, TYPES)
    ## use linestack (and exit) if only one axis was chosen, and
    ## display includes row or column scores
    if (length(choices) == 1) {
        ## Only one set of scores: plot them
        if (length(g) == 1)
            pl <- linestack(g[[1]], ...)
        ## The order of scores is species, sites, constraints, biplot,
        ## centroids: plot two first in this order, but species scores
        ## on the left
        else {
            hasSpec <- names(g)[1] == "species"
            ylim <- range(c(g[[1]], g[[2]]), na.rm = TRUE)
            pl <- linestack(g[[1]], ylim = ylim,
                            side = ifelse(hasSpec, "left", "right"), ...)
            linestack(g[[2]], ylim = ylim,
                      side = ifelse(hasSpec, "right","left"),
                      add = TRUE, ...)
            }
        return(invisible(pl))
    }
    if (missing(xlim)) {
        xlim <- range(0, g$species[, 1], g$sites[, 1], g$constraints[, 1],
                      g$biplot[, 1], g$regression[,1],
                      if (length(g$centroids) > 0 && all(is.na(g$centroids))) NA else g$centroids[, 1],
                      na.rm = TRUE)
    }
    if (!any(is.finite(xlim)))
        stop("no finite scores to plot")
    if (missing(ylim)) {
        ylim <- range(0, g$species[, 2], g$sites[, 2], g$constraints[, 2],
                      g$biplot[, 2], g$regression[,2],
                      if (length(g$centroids) > 0 && all(is.na(g$centroids))) NA else g$centroids[, 2],
                      na.rm = TRUE)
    }
    plot(g[[1]], xlim = xlim, ylim = ylim, type = "n", asp = 1,
         ...)
    abline(h = 0, lty = 3)
    abline(v = 0, lty = 3)
    if (!is.null(g$species)) {
        if (type == "text")
            text(g$species, rownames(g$species), col = "red",
                 cex = cex)
        else if (type == "points")
            points(g$species, pch = "+", col = "red", cex = cex)
    }
    if (!is.null(g$sites)) {
        if (type == "text")
            text(g$sites, rownames(g$sites), cex = cex)
        else if (type == "points")
            points(g$sites, pch = 1, cex = cex)
    }
    if (!is.null(g$constraints)) {
        if (type == "text")
            text(g$constraints, rownames(g$constraints), cex = cex,
                 col = "darkgreen")
        else if (type == "points")
            points(g$constraints, pch = 2, cex = cex, col = "darkgreen")
    }
    if (!is.null(g$biplot) && nrow(g$biplot) > 0 && type != "none") {
        if (length(display) > 1) {
            mul <- ordiArrowMul(g$biplot)
        }
        else mul <- 1
        attr(g$biplot, "arrow.mul") <- mul
        arrows(0, 0, mul * g$biplot[, 1], mul * g$biplot[, 2],
               length = 0.05, col = "blue")
        biplabs <- ordiArrowTextXY(mul * g$biplot, rownames(g$biplot),
                                   cex = cex, rescale = FALSE)
        text(biplabs, rownames(g$biplot), col = "blue", cex = cex)
    }
    if (!is.null(g$regression) && nrow(g$regression > 0) && type != "none") {
        rcol <- "purple4"
        if (length(display) > 1) {
            mul <- ordiArrowMul(g$regression)
        }
        else mul <- 1
        attr(g$regression, "arrow.mul") <- mul
        arrows(0, 0, mul * g$regression[, 1], mul * g$regression[, 2],
               length = 0.05, col = rcol)
        biplabs <- ordiArrowTextXY(mul * g$regression, rownames(g$regression),
                                   cex = cex, rescale = FALSE)
        text(biplabs, rownames(g$regression), col = rcol, cex = cex)
    }
    if (!is.null(g$centroids) && !anyNA(g$centroids) && type !=
        "none") {
        if (type == "text")
            text(g$centroids, rownames(g$centroids), col = "blue", cex = cex)
        else if (type == "points")
            points(g$centroids, pch = "x", col = "blue", cex = cex)
    }
    if (!is.null(g$default) && type != "none") {
        if (type == "text")
            text(g$default, rownames(g$default), cex = cex)
        else if (type == "points")
            points(g$default, pch = 1, cex = cex)
    }
    class(g) <- "ordiplot"
    invisible(g)
}

## vegan::plot.rda needed because klaR::plot.rda would be used
## instead if klaR package is loaded

`plot.rda`<-
    function(x, ...)
{
    ## not vegan rda?
    if (!("CA" %in% names(x)))
        stop(gettextf("%s is not a vegan rda object",
                      sQuote(deparse(substitute(x)))))
    NextMethod("plot", x, ...)
}
