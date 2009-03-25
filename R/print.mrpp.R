"print.mrpp" <-
function (x, digits = max(3, getOption("digits") - 3), ...) 
{
### A print function for mrpp objects
### x -- An object of class "mrpp."
#### cat = print
    cat("\nCall:\n")
    cat(deparse(x$call), "\n\n") 
    cat("Dissimilarity index:", x$distance, "\n")
    cat("Weights for groups: ", switch(x$weight.type, "n", "n-1", "n(n-1)", "n(n-1)/2"), "\n\n")
    cat("Class means and counts:\n\n")
    print(noquote(rbind("delta" = formatC(x$classdelta, digits = digits),
                        "n" = formatC(x$n, digits=0))))
    cat("\n")
    if (!is.na(x$CS)) {
        cat("Classification strength: ")
        cat(formatC(x$CS, digits = digits), "\n")
    }
    cat("Chance corrected within-group agreement A: ")
    cat(formatC(x$A, digits = digits), "\n")
    cat("Based on observed delta", formatC(x$delta), "and expected delta",
        formatC(x$E.delta),"\n\n")
    nperm <- x$permutations
    if (nperm) {
        cat("Significance of delta:", format.pval(x$Pvalue, eps = 1/nperm), 
            "\n")
        cat("Based on ", nperm, " permutations")
    }
    if (!is.null(x$strata)) 
        cat(", stratified within", x$strata)
    cat("\n\n")
    invisible(x)
}
