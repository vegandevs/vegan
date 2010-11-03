`plot.prc` <-
    function (x, species = TRUE, select, scaling = 3, axis = 1, type = "l",
              xlab, ylab, ylim, lty = 1:5, col = 1:6, pch, legpos, cex = 0.8,
              ...)
{
    ## save level names before getting the summary
    levs <- x$terminfo$xlev[[2]]
    x <- summary(x, scaling = scaling, axis = axis)
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    b <- t(coef(x))
    xax <- rownames(b)
    if (missing(xlab)) 
        xlab <- x$names[1]
    if (missing(ylab))
        ylab <- "Effect"
    if (!missing(select))
        x$sp <- x$sp[select]
    if (missing(ylim))
        ylim <- if (species)
            range(b, x$sp, na.rm = TRUE)
        else range(b, na.rm = TRUE)
    if (species) {
        op <- par("mai")
        mrg <- max(strwidth(names(x$sp), cex = cex, units = "in")) +
            strwidth("mmm", cex = cex, units = "in")
        par(mai = c(op[1:3], max(op[4], mrg)))
    }
    if (missing(pch))
        pch <- as.character(1:nrow(b))
    matplot(xax, b, type = type, xlab = xlab, ylab = ylab, ylim = ylim,
            cex = cex, lty = lty, col = col, pch = pch, ...)
    abline(h = 0, col = "gray")
    if (species) {
        linestack(x$sp, at = par("usr")[2], add = TRUE, hoff = 1,
                  cex = cex, ...)
        rug(x$sp, side = 4)
    }
    if (missing(legpos)) {
        holes <- abs(par("usr")[3:4] - range(b, na.rm = TRUE))
        if (holes[1] > holes[2])
            legpos <- "bottomleft"
        else legpos <- "topleft"
    }
    if (!is.na(legpos)) {
        nl <- length(levs)
        pp <- type %in% c("b", "p")
        pl <- type %in% c("b", "l")
        if (length(lty) == 1)
            lty <- rep(lty, nl-1)
        legend(legpos, legend = levs, col = c("gray", col),
               lty = if (pl) lty[c(1,1:(nl-1))],
               pch = if (pp) pch, cex = cex, title = x$names[2])
    }
    invisible()
}
