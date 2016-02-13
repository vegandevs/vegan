"print.radline" <-
    function (x, ...) 
{
    cat("\nRAD model:", x$model, "\n")
    cat("Family:", family(x)$family, "\n")
    cat("No. of species: ", length(x$y), "\nTotal abundance:", 
        sum(x$y), "\n\n")
    p <- coef(x)
    dev <- deviance(x)
    AIC <- AIC(x)
    BIC <- AIC(x, k = log(length(x$y)))
    tmp <- c(p, dev, AIC, BIC)
    names(tmp) <- c(names(p), "Deviance", "AIC", "BIC")
    print(tmp, ...)
    invisible(x)
}
