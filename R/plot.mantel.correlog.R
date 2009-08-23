`plot.mantel.correlog` <-
    function(x, alpha=0.05, ...)
{
    lim <- max(x$n.tests)
    plot(x$mantel.res[1:lim,1],x$mantel.res[1:lim,3],
         xlab="Distance class index", ylab="Mantel correlation", pch=22)
    if(x$mult == "none") {
	signif <- which((x$mantel.res[1:lim,4] <= alpha))
    } else {
	signif <- which((x$mantel.res[1:lim,5] <= alpha))
    }
    lines(x$mantel.res[1:lim,1], x$mantel.res[1:lim,3])
    points(x$mantel.res[1:lim,1], x$mantel.res[1:lim,3], pch=22, bg="white")
    points(x$mantel.res[signif,1], x$mantel.res[signif,3], pch=22, bg="black")
    abline(a=0, b=0, col="red")
    invisible()
}
