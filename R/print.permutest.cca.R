`print.permutest.cca` <-
    function (x, ...)
{
    EPS <- sqrt(.Machine$double.eps)
    cat("\nPermutation test for", x$method, "under", x$model, "model", "\n\n")
    if (x$nperm > 0)
        cat(howHead(x$control), "\n")
    writeLines(strwrap(pasteCall(x$testcall, prefix = "Model:")))
    if (x$nperm > 0)
        Pval <- (colSums(sweep(x$F.perm, 2, x$F.0 - EPS, ">=")) + 1)/(x$nperm + 1)
    else
        Pval <- NA
    cat("Permutation test for ")
    if (x$first)
        cat("first constrained eigenvalue\n")
    else if (length(x$df) <= 2)
        cat("all constrained eigenvalues\n")
    else
        cat("sequential contrasts\n")
    anotab <- data.frame(x$df, x$chi, c(x$F.0, NA), c(Pval, NA))
    colnames(anotab) <- c("Df", "Inertia", "F", "Pr(>F)")
    rownames(anotab) <- c(x$termlabels, "Residual")
    class(anotab) <- c("anova", "data.frame")
    print(anotab)
    invisible(x)
}
