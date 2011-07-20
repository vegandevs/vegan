`ordiplot3d` <-
    function (object, display = "sites", choices = 1:3, ax.col = 2, 
              arr.len = 0.1, arr.col = 4, envfit, xlab, ylab, zlab, ...) 
{
    require(scatterplot3d) || stop("Requires package 'scatterplot3d'")
    x <- scores(object, display = display, choices = choices, ...)
    if (missing(xlab)) xlab <- colnames(x)[1]
    if (missing(ylab)) ylab <- colnames(x)[2]
    if (missing(zlab)) zlab <- colnames(x)[3]
    pl <- ordiArgAbsorber(x[, 1], x[, 2], x[, 3],  
                        xlab = xlab, ylab = ylab, zlab = zlab,
                          FUN = "scatterplot3d", ...)
    pl$points3d(range(x[, 1]), c(0, 0), c(0, 0), type = "l", 
                col = ax.col)
    pl$points3d(c(0, 0), range(x[, 2]), c(0, 0), type = "l", 
                col = ax.col)
    pl$points3d(c(0, 0), c(0, 0), range(x[, 3]), type = "l", 
                col = ax.col)
    if (!missing(envfit) ||
        (!is.null(object$CCA) && object$CCA$rank > 0)) {
        if (!missing(envfit)) 
            object <- envfit
        bp <- scores(object, dis = "bp", choices = choices, ...)
        cn <- scores(object, dis = "cn", choices = choices, ...)
        if (!is.null(cn) && !any(is.na(cn))) {
            bp <- bp[!(rownames(bp) %in% rownames(cn)), , drop = FALSE]
            cn.xyz <- pl$xyz.convert(cn)
            points(cn.xyz, pch = "+", cex = 2, col = arr.col)
        }
        if (!is.null(bp) && nrow(bp) > 0) {
            tmp <- pl$xyz.convert(bp)
            mul <- ordiArrowMul(cbind(tmp$x, tmp$y), fill=1)
            bp.xyz <- pl$xyz.convert(bp * mul)
            orig <- pl$xyz.convert(0, 0, 0)
            arrows(orig$x, orig$y, bp.xyz$x, bp.xyz$y, length = arr.len, 
                   col = arr.col)
        }
    }
    tmp <- pl$xyz.convert(x)
    pl$points <- cbind(tmp$x, tmp$y)
    rownames(pl$points) <- rownames(x)
    if (exists("bp.xyz")) {
        pl$arrows <- cbind(bp.xyz$x, bp.xyz$y)
        rownames(pl$arrows) <- rownames(bp)
    }
    if (exists("cn.xyz")) {
        pl$centroids <- cbind(cn.xyz$x, cn.xyz$y)
        rownames(pl$centroids) <- rownames(cn)
    }
    class(pl) <- c("ordiplot3d", "ordiplot")
    invisible(pl)
}
