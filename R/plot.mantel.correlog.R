'plot.mantel.correlog' <- function(x, ...)
{
lim = max(x$n.tests)
plot(x$mantel.res[1:lim,1],x$mantel.res[1:lim,3], xlab="Distance class index", ylab="Mantel r", pch=22)
if(x$mult=="none") {
	signif = which((x$mantel.res[1:lim,4] <= 0.05))
	} else {
	signif = which((x$mantel.res[1:lim,5] <= 0.05))
	}
points(x$mantel.res[signif,1], x$mantel.res[signif,3], pch=22, bg="black")
abline(a=0, b=0, col="red") 
invisible()
}