"plot.fisherfit" <-
    function(x, xlab = "Frequency", ylab = "Species", bar.col = "skyblue",
             line.col= "red", lwd=2, ...)
{
    freq <- as.numeric(names(x$fisher))
    plot(freq, x$fisher, ylab=ylab, xlab=xlab,
         ylim=c(0,max(x$fisher)),  xlim=c(0.5, max(freq)+0.5), type="n", ...)
    rect(freq-0.5,0,freq+0.5,x$fisher, col=bar.col, ...)
    alpha <- x$estimate
    k <- x$nuisance
    curve(alpha*k^x/x, 1, max(freq), col=line.col, lwd=lwd, add=TRUE, ...)
    invisible()
}
