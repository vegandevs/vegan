`print.isomap` <-
function (x, ...) 
{
	cat("\nIsometric Feature Mapping (isomap)\n\n")
	cat("Call:\n")
	cat(deparse(x$call), "\n\n")
	cat("Distance method:", x$method, "\n")
	cat("Criterion:", x$criterion, "=", x$critval, "\n")
	if(!is.null(x$take))
	    cat("Data were fragmented, analysed", sum(x$take), "of", length(x$take), "points\n")
	invisible(x)
}

