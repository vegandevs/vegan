`summary.isomap` <-
function (object, ...)
{
	axes <- min(axes, ncol(object$points))
	out <- list()
	out$call <- object$call
	out$net <- object$net
	n <- nrow(object$points)
	out$ndis <- n * (n-1) / 2
	out$nnet <- nrow(object$net)
	class(out) <- "summary.isomap"
	out
}

