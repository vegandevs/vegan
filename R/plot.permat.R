## S3 plot method for permat
`plot.permat` <-
function(x, type = "bray", ylab, xlab, col, lty, lowess=TRUE, plot=TRUE, text=TRUE, ...)
{
    type <- match.arg(type, c("bray", "chisq"))
    if (missing(xlab))
        xlab <- "Runs"
    if (missing(col))
        col <- c(2,4)
    if (missing(lty))
        lty <- c(1,2)
    n <- attr(x, "times")
    toplot <- numeric(n)
    if (type == "bray") {
        toplot <- summary(x)$bray
        if (missing(ylab))
            ylab <- "Bray-Curtis dissimilarity"
    }
    if (type == "chisq") {
        toplot <- summary(x)$chisq
        if (missing(ylab))
            ylab <- "Chi-squared"
    }
    if (plot) {
        plot(toplot,type="n",ylab=ylab,xlab=xlab, ...)
        lines(toplot,col=col[1], lty=lty[1])
        if (lowess)
            lines(lowess(toplot),col=col[2], lty=lty[2])
        if (text) title(sub=paste("(mean = ", substitute(z, list(z=round(mean(toplot),3))), 
            ", min = ", substitute(z, list(z=round(min(toplot),3))),
            ", max = ", substitute(z, list(z=round(max(toplot),3))), ")", sep=""))
    }
    invisible(toplot)
}
