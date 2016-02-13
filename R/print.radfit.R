"print.radfit" <-
function(x,  digits = max(3, getOption("digits") - 2),  ...) 
{
    cat("\nRAD models, family", x$family$family, "\n")
    cat("No. of species ", length(x$y), ", total abundance ", 
        sum(x$y), "\n\n", sep = "")
    p <- coef(x)
    if (any(!is.na(p)))
        p <- formatC(p, format="g", flag = " ", digits = digits)
    p <- apply(p, 2, function(x) gsub("NA", " ", x))
    aic <- sapply(x$models, AIC)
    bic <- sapply(x$models, AIC, k = log(length(x$y)))
    dev <- sapply(x$models, deviance)
    stats <- format(cbind(Deviance = dev, AIC = aic, BIC = bic), digits = digits, ...)
    out <- cbind(p, stats)
    print(out, quote=FALSE)
    invisible(x)
}
