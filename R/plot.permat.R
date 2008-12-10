## S3 plot method for permat
`plot.permat` <-
function(x, ylab, xlab, col, lty, plot=TRUE, text=TRUE, ...)
{
    if (missing(ylab))
        ylab <- "Bray-Curtis dissimilarity"
    if (missing(xlab))
        xlab <- "Runs"
    if (missing(col))
        col <- c(2,4)
    if (missing(lty))
        lty <- c(1,2)
    n <- attr(x, "times")
    bray <- numeric(n)
    bray <- summary(x)$bray
    if (plot) {
        plot(bray,type="n",ylab=ylab,xlab=xlab, ...)
        lines(bray,col=col[1], lty=lty[1])
        lines(lowess(bray),col=col[2], lty=lty[2])
        if (text) title(sub=paste("(mean = ", substitute(z, list(z=round(mean(bray),3))), 
            ", min = ", substitute(z, list(z=round(min(bray),3))),
            ", max = ", substitute(z, list(z=round(max(bray),3))), ")", sep=""))
    }
    invisible(bray)
}
