"print.summary.procrustes" <-
function (x, digits = x$digits, ...) 
{
    cat("\nCall:\n")
    cat(deparse(x$call), "\n")
    cat("\nNumber of objects:", x$n, "   Number of dimensions:", 
        x$k, "\n")
    cat("\nProcrustes sum of squares:  ")
    cat("\n", formatC(x$ss, digits = digits), "\n")
    cat("Procrustes root mean squared error: ")
    cat("\n", formatC(x$rmse, digits = digits), "\n")
    cat("Quantiles of Procrustes errors:\n")
    nam <- c("Min", "1Q", "Median", "3Q", "Max")
    rq <- structure(quantile(x$resid), names = nam)
    print(rq, digits = digits, ...)
    cat("\nRotation matrix:\n")
    print(x$rotation, digits = digits, ...)
    cat("\nTranslation of averages:\n")
    print(x$translation, digits = digits, ...)
    cat("\nScaling of target:\n")
    print(x$scale, digits = digits, ...)
    cat("\n")
    invisible(x)
}

