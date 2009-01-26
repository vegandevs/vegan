"ordirgl" <-
    function (object, display = "sites", choices = 1:3, type = "p", 
              ax.col = "red", arr.col = "yellow", text, envfit, ...) 
{
    if (!require(rgl)) 
        stop("Requires package 'rgl'")
    oldpak <- compareVersion(packageDescription("rgl", field="Version"), "0.65") == -1
    x <- scores(object, display = display, choices = choices, 
                ...)
    if (ncol(x) < 3) 
        stop("3D display needs three dimensions...")
    rgl.clear()
    if (type == "p") 
        rgl.points(x[, 1], x[, 2], x[, 3], ...)
    else if (type == "t") {
        if (missing(text)) 
            text <- rownames(x)
        rgl.texts(x[, 1], x[, 2], x[, 3], text, ...,
                  if (oldpak) justify = "center" else adj = 0.5)
    }
    rgl.lines(range(x[, 1]), c(0, 0), c(0, 0), col = ax.col)
    rgl.lines(c(0, 0), range(x[, 2]), c(0, 0), col = ax.col)
    rgl.lines(c(0, 0), c(0, 0), range(x[, 3]), col = ax.col)
    rgl.texts(1.1 * max(x[, 1]), 0, 0, colnames(x)[1], col = ax.col, 
              if (oldpak) justify = "center" else adj = 0.5)
    rgl.texts(0, 1.1 * max(x[, 2]), 0, colnames(x)[2], col = ax.col, 
              if (oldpak) justify = "center" else adj = 0.5)
    rgl.texts(0, 0, 1.1 * max(x[, 3]), colnames(x)[3], col = ax.col, 
              if (oldpak) justify = "center" else adj = 0.5)
    if (!missing(envfit) || !is.null(object$CCA)) {
        if (!missing(envfit)) 
            object <- envfit
        bp <- scores(object, dis = "bp", choices = choices)
        cn <- scores(object, dis = "cn", choices = choices)
        if (!is.null(cn) && !any(is.na(cn))) {
            bp <- bp[!(rownames(bp) %in% rownames(cn)), , drop = FALSE]
            rgl.texts(cn[, 1], cn[, 2], cn[, 3], rownames(cn), 
                      col = arr.col, if (oldpak) justify = "center" else adj = 0.5)
            rgl.points(cn[, 1], cn[, 2], cn[, 3], size = 5, col = arr.col)
        }
        if (!is.null(bp) && nrow(bp) > 0) {
            mul <- c(range(x[, 1]), range(x[, 2]), range(x[, 
                                                           3]))/c(range(bp[, 1]), range(bp[, 2]), range(bp[, 
                                                                                                           3]))
            mul <- mul[is.finite(mul) & mul > 0]
            mul <- min(mul)
            bp <- bp * mul
            for (i in 1:nrow(bp)) {
                rgl.lines(c(0, bp[i, 1]), c(0, bp[i, 2]), c(0, 
                                                            bp[i, 3]), col = arr.col)
                rgl.texts(1.1 * bp[i, 1], 1.1 * bp[i, 2], 1.1 * 
                          bp[i, 3], rownames(bp)[i], col = arr.col,
                          if (oldpak) justify = "center" else adj = 0.5)
            }
        }
    }
    invisible()
}
