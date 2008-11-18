"plot.radfit" <-
    function (x, BIC = FALSE, legend = TRUE, ...) 
{
    if (length(x$y) == 0)
        stop("No species, nothing to plot")
    out <- plot(x$y, ...)
    if (length(x$y) == 1)
        return(invisible(out))
    fv <- fitted(x)
    if (BIC) 
        k = log(length(x$y))
    else k = 2
    emph <- which.min(sapply(x$models, AIC, k = k))
    lwd <- rep(1, ncol(fv))
    lwd[emph] <- 3
    matlines(fv, lty = 1, lwd = lwd, ...)
    if (legend) {
        nm <- names(x$models)
        legend("topright", legend = nm, lty = 1, lwd = lwd, col = 1:6)
    }
    invisible(out)
}
