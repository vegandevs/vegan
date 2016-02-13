`print.varpart234` <-
    function(x, digits = 5, ...)
{
    cat("No. of explanatory tables:", x$nsets, "\n")
    cat("Total variation (SS):", format(x$SS.Y, digits=digits), "\n")
    if (x$ordination == "rda")
        cat("            Variance:", format(x$SS.Y/(x$n-1), digits=digits), "\n")
    cat("No. of observations:",  x$n, "\n")
    cat("\nPartition table:\n")
    out <- rbind(x$fract, "Individual fractions" = NA, x$indfract)
    if (x$nsets > 3)
        out <- rbind(out, "Controlling 2 tables X" = NA, x$contr2)
    if (x$nsets > 2)
        out <- rbind(out, "Controlling 1 table X" = NA, x$contr1)
    out[,2:3] <- round(out[,2:3], digits=digits)
    out[,1:4] <- sapply(out[,1:4], function(x) gsub("NA", "  ", format(x, digits=digits)))
    print(out)
    cat("---\nUse function", sQuote(x$ordination),
        "to test significance of fractions of interest\n")
    if (!is.null(x$bigwarning))
        for (i in seq_along(x$bigwarning))
            warning("collinearity detected: redundant variable(s)  between tables ",
                    x$bigwarning[i],
                    "\nresults are probably incorrect: remove redundant variable(s) and repeat the analysis",
                    call. = FALSE)
    invisible(x)
}

