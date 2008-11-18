`plot.isomap` <-
function (x, net=TRUE, n.col = "gray", ...) 
{
	pl <- ordiplot(x, display="sites", ...)
	if (net) {
		z <- scores(pl, "sites")
		k <- x$net
		segments(z[k[,1],1], z[k[,1],2], z[k[,2],1], z[k[,2],2], col=n.col)
	}
	invisible(pl)
}

