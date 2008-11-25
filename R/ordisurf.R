"ordisurf" <-
    function (x, y, choices = c(1, 2), knots = 10, family = "gaussian", 
              col = "red", thinplate = TRUE, add = FALSE, display = "sites", 
              w = weights(x), main, nlevels = 10, levels, labcex = 0.6,
              bubble = FALSE, cex = 1, ...) 
{
    weights.default <- function(object, ...) NULL
    GRID = 25
    w <- eval(w)
    if (!is.null(w) && length(w) == 1) 
        w <- NULL
    if (!require(mgcv)) 
        stop("Requires package `mgcv'")
    X <- scores(x, choices = choices, display = display, ...)
    x1 <- X[, 1]
    x2 <- X[, 2]
    if (thinplate) 
        mod <- gam(y ~ s(x1, x2, k = knots), family = family, 
                   weights = w)
    else mod <- gam(y ~ s(x1, k = knots) + s(x2, k = knots), 
                    family = family, weights = w)
    xn1 <- seq(min(x1), max(x1), len=GRID)
    xn2 <- seq(min(x2), max(x2), len=GRID)
    newd <- expand.grid(x1 = xn1, x2 = xn2)
    fit <- predict(mod, type = "response", newdata=as.data.frame(newd))
    poly <- chull(cbind(x1,x2))
    poly <- c(poly, poly[1])
    npol <- length(poly)
    np <- nrow(newd)
    inpoly <- numeric(np)
    inpoly <- .C("pnpoly", as.integer(npol), as.double(x1[poly]), as.double(x2[poly]),
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
            main <- deparse(substitute(y))
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
