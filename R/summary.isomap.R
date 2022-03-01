`summary.isomap` <-
function (object, axes=4, ...) 
{
	axes <- min(axes, ncol(object$points))
	out <- list()
	out$call <- object$call
	out$points <- object$points[,1:axes]
	out$net <- object$net
	n <- nrow(object$points)
	out$ndis <- n * (n-1) / 2
	out$nnet <- nrow(object$net)
	class(out) <- "summary.isomap"
	out
}

