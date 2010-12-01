"plot.procrustes" <-
    function (x, kind = 1, choices = c(1,2), to.target = TRUE,
              type = "p", xlab, ylab, main, ar.col = "blue",
              len = 0.05, cex = 0.7,  ...) 
{
    type <- match.arg(type, c("points", "text", "none"))
    if (to.target) {
        tails <- x$Yrot[, choices]
        heads <- x$X[, choices]
    } else {
        tails <- x$X[, choices]
        heads <- x$Yrot[, choices]
    }
    if (missing(main)) 
        main <- "Procrustes errors"
    if (kind <= 1) {
        formals(arrows) <- c(formals(arrows), alist(... = ))
        if (missing(xlab)) 
            xlab <- paste("Dimension", choices[1])
        if (missing(ylab)) 
            ylab <- paste("Dimension", choices[2])
        xrange <- range(tails[, 1], heads[, 1])
        yrange <- range(tails[, 2], heads[, 2])
        plot(xrange, yrange, xlab = xlab, ylab = ylab, main = main, 
             type = "n", asp = 1, ...)
        if (kind > 0) {
            abline(v = 0, lty = 2)
            abline(h = 0, lty = 2)
            if (ncol(x$rotation) == 2) {
                ## Draw rotated axes only if they visibly differ from
                ## unrotated axes
                b <- abs(x$rotation[1,1])
                EPS <- 1e-8
                if (1 - b > EPS && b > EPS) {
                    abline(0, tan(acos(x$rotation[1, 1])), lty = 1)
                    abline(0, 1/tan(acos(-x$rotation[1, 1])), lty = 1)
                } else {
                    abline(v = 0, h = 0)
                }
            }
            else {
            	Y <- x$Yrot %*% t(x$rotation)
            	for (k in 1:ncol(Y)) {
                    tmp <- matrix(0, nrow=2, ncol = ncol(Y))
                    tmp[,k] <- range(Y[,k])
                    tmp <- tmp %*% x$rotation
                    lines(tmp[,choices], lty=1)
                    text(tmp[2,choices[1]], tmp[2,choices[2]], as.character(k))
                }
            }
            if (type != "none") {
                ow <- options(warn = -1)
                arrows(tails[, 1], tails[, 2], heads[, 1], heads[, 2],
                       col = ar.col, len = len, ...)
                options(ow)
                if (type == "text" && !is.null(rownames(tails)))
                    ordilabel(tails, cex = cex, ...)
                else
                    points(tails, cex = cex, ...)
            }
        }
        out <- list(heads = heads, points = tails)
        class(out) <- "ordiplot"
    }
    else if (kind == 2) {
        if (missing(xlab)) 
            xlab <- "Index"
        if (missing(ylab)) 
            ylab <- "Procrustes residual"
        res <- residuals(x)
        q <- quantile(res)
        plot(res, type = "h", xlab = xlab, ylab = ylab, main = main, 
             ...)
        abline(h = q[2:4], lty = c(2, 1, 2))
        out <- list(sites = cbind(seq(along = res), res))
        class(out) <- "ordiplot"
    }
    invisible(out)
}
