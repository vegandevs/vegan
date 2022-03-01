`print.summary.isomap` <-
function (x, ...) 
{
	cat("\nCall:\n")
	cat(deparse(x$call), "\n\n")
	cat("Points:\n")
	print(x$points, ...)
	cat("\nRetained dissimilarities between points:\n")
	print(t(x$net), ...)
	cat("\nRetained", x$nnet, "of", x$ndis, "dissimilarities\n")
	invisible(x)
}

