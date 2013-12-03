`print.permutest.cca` <-
    function (x, ...) 
{
    cat("\nPermutation test for", x$method, "\n\n")
    cat(howHead(x$control), "\n")
    writeLines(strwrap(pasteCall(x$testcall)))
    Pval <- (sum(x$F.perm >= x$F.0) + 1)/(x$nperm + 1)
    cat("Permutation test for ")
    if (x$first)
        cat("first constrained eigenvalue\n")
    else
        cat("all constrained eigenvalues\n")
    cat("Pseudo-F:\t", x$F.0, "(with", paste(x$df, collapse = ", "),
        "Degrees of Freedom)\n")
    cat("Significance:\t", format.pval(Pval), 
        "\n\n")
    invisible(x)
}
