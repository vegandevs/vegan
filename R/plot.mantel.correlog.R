'plot.mantel.correlog' <- function(x, ...)
{
how.many.classes = which(x$mantel.res[,3] != "NA")
lim = max(how.many.classes)
plot(x$mantel.res[1:lim,1],x$mantel.res[1:lim,3], xlab="Distance class index", ylab="Mantel r", pch=22)
not.signif = which((x$mantel.res[1:lim,4] <= 0.05))
points(x$mantel.res[not.signif,1], x$mantel.res[not.signif,3], pch=22, bg="black")
abline(a=0, b=0, col="red") 
invisible()
}
