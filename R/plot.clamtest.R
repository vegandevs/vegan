plot.clamtest <- function(x, xlab, ylab, main,
    pch=21:24, col.points=1:4, col.lines=2:4, lty=1:3,
    position="bottomright", ...) {
    summ <- summary(x)
    glabel <- summ$labels
    if (missing(main))
        main <- "Species Classification"
    if (missing(xlab))
        xlab <- paste(glabel[2], "(abundance + 1)")
    if (missing(ylab))
        ylab <- paste(glabel[1], "(abundance + 1)")
    Y <- x[,2]
    X <- x[,3]
    minval <- summ$minv
    ## plot the dots
    rr <- range(X+1,Y+1)
    plot(X+1, Y+1, log = "xy", xaxt = "n", yaxt = "n",
        col=col.points[as.integer(x$Classes)],
        pch=pch[as.integer(x$Classes)], 
        xlab=xlab, ylab=ylab, main=main,
        xlim=rr, ylim=rr, ...)
    axis(1, c(1,10,100,1000,10000))
    axis(2, c(1,10,100,1000,10000))
    ## too rare threshold
    Ymin <- minval[[1]][1,2]
    Xmin <- minval[[2]][1,1]
    lines(rep(Xmin, 2)+1, c(0, 1)+1, col=col.lines[1], lty=lty[1])
    lines(c(0, 1)+1, rep(Ymin, 2)+1, col=col.lines[1], lty=lty[1])
    tmp <- approx(c(Xmin, 1), c(1, Ymin))
    lines(tmp$x+1, tmp$y+1, col=col.lines[1], lty=lty[1])
    ## Y vs. gen threshold
    lines(minval[[1]]+1, col=col.lines[2], lty=lty[2])
    ## X vs. gen threshold
    lines(minval[[2]]+1, col=col.lines[3], lty=lty[3])
    if (!is.null(position))
        legend(position, col=col.points, pch=pch, 
            legend=rownames(summ$summary))
    invisible(x)
}
