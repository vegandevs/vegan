`print.summary.isomap` <-
function (x, ...)
{
	cat("\nCall:\n")
	cat(deparse(x$call), "\n\n")
	cat("\nRetained dissimilarities between points:\n")
	prmatrix(t(x$net), collab = rep("", x$nnet) , ...)
	cat("\nRetained", x$nnet, "of", x$ndis, "dissimilarities\n")
	invisible(x)
}

