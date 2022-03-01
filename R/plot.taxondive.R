`plot.taxondive` <-
function (x, ...) 
{
	plot(x$Species, x$Dplus, xlab="Number of Species", ylab=expression(Delta^"+"), ...) 
	i <- order(x$Species)
	abline(h=x$EDplus, ...)
	lines(x$Species[i], x$EDplus - 2*x$sd.Dplus[i], ...)
     lines(x$Species[i], x$EDplus + 2*x$sd.Dplus[i], ...)
}

