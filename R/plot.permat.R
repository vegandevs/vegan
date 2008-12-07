## S3 plot method for permat
`plot.permat` <-
function(x, ylab="Bray-Curtis dissimilarity", xlab="Runs", col=c(2,4), lty=c(1,2), plot=TRUE, text=TRUE, ...)
{
    n <- attr(x, "times")
    bray <- numeric(n)
#    for (i in 1:n) bray[i] <- sum(abs(x$orig-x$perm[[i]]))/sum(x$orig+x$perm[[i]])
    bray <- sapply(x$perm, function(z) sum(abs(x$orig - z)) / sum(x$orig + z))
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
