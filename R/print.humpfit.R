"print.humpfit" <-
    function(x, ...)
{
    cat("\nHump-backed Null model of richness vs. productivity\n\n")
    cat("Family:", family(x)$family,"\n")
    cat("Link function: Fisher diversity\n\n")
    cat("Coefficients:\n\n")
    print(coef(x))
    cat("\nDeviance", deviance(x), "with", df.residual(x))
    cat(" residual degrees of freedom\n")
    invisible(x)
}
