"ordisurf" <-
    function (x, y, choices = c(1, 2), knots = 10, family = "gaussian", 
              col = "red", thinplate = TRUE, add = FALSE, display = "sites", 
              w = weights(x), main, nlevels = 10, levels, labcex = 0.6,
              bubble = FALSE, cex = 1, ...) 
{
    weights.default <- function(object, ...) NULL
    GRID = 31
    w <- eval(w)
    if (!is.null(w) && length(w) == 1) 
        w <- NULL
    require(mgcv)  || stop("Requires package 'mgcv'")
    X <- scores(x, choices = choices, display = display, ...)
    ## The original name of 'y' may be lost in handling NA: save for
    ## plots
    yname <- deparse(substitute(y))
    kk <- complete.cases(X) & !is.na(y)
    if (!all(kk)) {
        X <- X[kk, , drop = FALSE]
        y <- y[kk]
        w <- w[kk]
    }
    x1 <- X[, 1]
    x2 <- X[, 2]
    if (knots <= 0)
        mod <- gam(y ~ x1 + x2, family = family, weights = w)
    else if (knots == 1)
        mod <- gam(y ~ poly(x1, 1) + poly(x2, 1),
                   family = family, weights = w)
    else if (knots == 2)
        mod <- gam(y ~ poly(x1, 2) + poly(x2, 2) + poly(x1, 1):poly(x2, 1),
                   family = family, weights = w)
    else if (thinplate) 
        mod <- gam(y ~ s(x1, x2, k = knots), family = family, 
                   weights = w)
    else mod <- gam(y ~ s(x1, k = knots) + s(x2, k = knots), 
                    family = family, weights = w)
    xn1 <- seq(min(x1), max(x1), len=GRID)
    xn2 <- seq(min(x2), max(x2), len=GRID)
    newd <- expand.grid(x1 = xn1, x2 = xn2)
    fit <- predict(mod, type = "response", newdata=as.data.frame(newd))
    poly <- chull(cbind(x1,x2))
    ## Move out points of the convex hull to have contour for all data
    ## points
    xhull1 <- x1[poly] + sign(x1[poly] - mean(x1[poly])) * diff(range(x1))/(GRID - 1)
    xhull2 <- x2[poly] + sign(x2[poly] - mean(x2[poly])) * diff(range(x2))/(GRID - 1)
    npol <- length(poly)
    np <- nrow(newd)
    inpoly <- numeric(np)
    inpoly <- .C("pnpoly", as.integer(npol), as.double(xhull1), as.double(xhull2),
                 as.integer(np), as.double(newd[,1]), as.double(newd[,2]),
                 inpoly = as.integer(inpoly), PACKAGE="vegan")$inpoly
    is.na(fit) <- inpoly == 0
    if (!add) {
        if (bubble) {
            if (is.numeric(bubble))
                cex <- bubble
            cex <- (y -  min(y))/diff(range(y)) * (cex-0.4) + 0.4
        }
        plot(X, asp = 1, cex = cex, ...)
    }
    if (!missing(main) || (missing(main) && !add)) {
        if (missing(main)) 
            main <- yname
        title(main = main)
    }
    if (missing(levels)) 
        levels <- pretty(range(fit, finite = TRUE), nlevels)
    
    contour(xn1, xn2, matrix(fit, nrow=GRID), col = col, add = TRUE,
            levels = levels, labcex = labcex,
            drawlabels = !is.null(labcex) && labcex > 0)
    mod$grid <- list(x = xn1, y = xn2, z = matrix(fit, nrow = GRID))
    class(mod) <- c("ordisurf", class(mod))
    return(mod)
}
