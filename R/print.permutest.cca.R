`print.permutest.cca` <-
    function (x, ...)
{
    EPS <- sqrt(.Machine$double.eps)
    cat("\nPermutation test for", x$method, "\n\n")
    cat(howHead(x$control), "\n")
    writeLines(strwrap(pasteCall(x$testcall)))
    Pval <- (colSums(sweep(x$F.perm, 2, x$F.0 - EPS, ">=")) + 1)/(x$nperm + 1)
    cat("Permutation test for ")
    if (x$first)
        cat("first constrained eigenvalue\n")
    else if (length(x$df) <= 2)
        cat("all constrained eigenvalues\n")
    else
        cat("sequential contrasts\n")
    anotab <- data.frame(x$df, x$chi, c(x$F.0, NA), c(Pval, NA))
    colnames(anotab) <- c("Df", "Inertia", "F", "Pr(>F)")
    if (nrow(anotab) == 2)
        rownames(anotab) <- c("Model", "Residual")
    else
        rownames(anotab)[nrow(anotab)] <- "Residual"
    class(anotab) <- c("anova", "data.frame")
    print(anotab)
    invisible(x)
}
