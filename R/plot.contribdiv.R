plot.contribdiv <-
function(x, sub, xlab, ylab, ylim, col, ...) {
    y <- x[,c(1,3)]
    if (missing(ylab))
        ylab <- paste("Diversity components (", attr(x, "index"), ")", sep = "")
    if (missing(xlab))
        xlab <- "Sites"
    if (missing(sub))
        sub <- paste("Differentiation coefficient = ", round(attr(x, "diff.coef"),3), sep = "")
    if (missing(ylim))
        ylim <- c(0, max(y))
    if (missing(col))
        col <- c("lightgrey", "darkgrey")
    matplot(y, type = "n", sub=sub, xlab=xlab, ylab=ylab, axes = FALSE,
            bty = "n", ...)
    polygon(c(1,1:nrow(y),nrow(y)), c(0,y$gamma,0), col=col[1])
    polygon(c(1,1:nrow(y),nrow(y)), c(0,y$alpha,0), col=col[2])
    axis(side = 1)
    axis(side = 2)
    box()
    invisible(x)
}
