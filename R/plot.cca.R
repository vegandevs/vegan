`plot.cca` <-
    function (x, choices = c(1, 2), display = c("sp", "wa", "cn"), 
              scaling = 2, type, xlim, ylim,  const, ...) 
{
    TYPES <- c("text", "points", "none")
    g <- scores(x, choices, display, scaling, const)
    if (!is.list(g)) 
        g <- list(default = g)
    ## Take care that there are names
    for (i in 1:length(g)) {
        if (length(dim(g[[i]])) > 1)
            rownames(g[[i]]) <- rownames(g[[i]], do.NULL = FALSE,
                                         prefix = substr(names(g)[i], 1, 3))
    }
    if (!is.null(g$centroids)) {
        if (is.null(g$biplot)) 
            g$biplot <- scores(x, choices, "bp", scaling)
        if (!is.na(g$centroids)[1]) {
            bipnam <- rownames(g$biplot)
            cntnam <- rownames(g$centroids)
            g$biplot <- g$biplot[!(bipnam %in% cntnam), , drop = FALSE]
            if (nrow(g$biplot) == 0) 
                g$biplot <- NULL
        }
    }
    if (missing(type)) {
        nitlimit <- 80
        nit <- max(nrow(g$spe), nrow(g$sit), nrow(g$con), nrow(g$def))
        if (nit > nitlimit) 
            type <- "points"
        else type <- "text"
    }
    else type <- match.arg(type, TYPES)
    if (missing(xlim))
        xlim <- range(g$spe[, 1], g$sit[, 1], g$con[, 1], g$default[, 
                                                                    1])
    if (missing(ylim))
        ylim <- range(g$spe[, 2], g$sit[, 2], g$con[, 2], g$default[, 
                                                                    2])
    plot(g[[1]], xlim = xlim, ylim = ylim, type = "n", asp = 1, 
         ...)
    abline(h = 0, lty = 3)
    abline(v = 0, lty = 3)
    if (!is.null(g$species)) {
        if (type == "text") 
            text(g$species, rownames(g$species), col = "red", 
                 cex = 0.7)
        else if (type == "points") 
            points(g$species, pch = "+", col = "red", cex = 0.7)
    }
    if (!is.null(g$sites)) {
        if (type == "text") 
            text(g$sites, rownames(g$sites), cex = 0.7)
        else if (type == "points") 
            points(g$sites, pch = 1, cex = 0.7)
    }
    if (!is.null(g$constraints)) {
        if (type == "text") 
            text(g$constraints, rownames(g$constraints), cex = 0.7, 
                 col = "darkgreen")
        else if (type == "points") 
            points(g$constraints, pch = 2, cex = 0.7, col = "darkgreen")
    }
    if (!is.null(g$biplot) && type != "none") {
        if (length(display) > 1) {
            mul <- ordiArrowMul(g$biplot)
        }
        else mul <- 1
        attr(g$biplot, "arrow.mul") <- mul
        arrows(0, 0, mul * g$biplot[, 1], mul * g$biplot[, 2], 
               len = 0.05, col = "blue")
        text(1.1 * mul * g$biplot, rownames(g$biplot), col = "blue")
        axis(3, at = c(-mul, 0, mul), labels = rep("", 3), col = "blue")
        axis(4, at = c(-mul, 0, mul), labels = c(-1, 0, 1), col = "blue")
    }
    if (!is.null(g$centroids) && !is.na(g$centroids) && type != 
        "none") {
        if (type == "text") 
            text(g$centroids, rownames(g$centroids), col = "blue")
        else if (type == "points") 
            points(g$centroids, pch = "x", col = "blue")
    }
    if (!is.null(g$default) && type != "none") {
        if (type == "text") 
            text(g$default, rownames(g$default), cex = 0.7)
        else if (type == "points") 
            points(g$default, pch = 1, cex = 0.7)
    }
    class(g) <- "ordiplot"
    invisible(g)
}
