"print.prestonfit" <-
    function (x, ...) 
{
    cat("\nPreston lognormal model\n")
    cat("Method:", x$method,"\n")
    cat("No. of species:", sum(x$freq), "\n\n")
    print(x$coefficients, ...)
    cat("\nFrequencies by Octave\n")
    print(rbind(Observed = x$freq, Fitted = x$fitted), ...)
    cat("\n")
    invisible(x)
}
