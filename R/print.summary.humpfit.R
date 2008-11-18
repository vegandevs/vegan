"print.summary.humpfit" <-
    function (x, ...) 
{
    cat("\nHump-backed Null model of richness vs. productivity\n\n")
    cat("Family:", x$family, "\n")
    cat("Link function: Fisher diversity\n\n")
    cat("Coefficients:\n\n")
    printCoefmat(x$est, ...)
    cat("\nDispersion parameter for", x$family, "family taken to be", x$dispersion,"\n")
    cat("\nDeviance", x$deviance, "with", x$df.residual)
    cat(" residual degrees of freedom\n")
    cat("AIC:", x$aic, "   BIC:", x$bic, "\n")
    cat("\nCorrelation of Coefficients:\n")
    correl <- format(round(x$correlation, 2), nsmall = 2)
    correl[!lower.tri(correl)] <- ""
    print(correl[-1, -3], quote = FALSE)
    cat("\nDiagnostics from nlm:\n")
    cat("Number of iterations: ", x$iter, ", code: ", x$code, 
        "\n", sep = "")
    invisible(x)
}
