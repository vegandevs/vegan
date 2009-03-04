## S3 print method for summary.permat
`print.summary.permat` <-
function(x, digits=2, ...)
{
    n <- attr(x$x, "times")
    cat("Summary of object of class 'permat'\n\nCall: ")
    print(x$x$call)
    cat("\nMatrix type:", attr(x$x, "mtype"), "\nPermutation type:", attr(x$x, "ptype"))
    if (attr(x$x, "ptype") == "swap") {
        cat("\nMethod: ", attr(x$x, "method"), sep = "")
        if (attr(x$x, "method") != "quasiswap") {
            cat(", burnin: ", attr(x$x, "burnin"), sep = "")
            cat(", thin: ", attr(x$x, "thin"), sep = "")
        }
    }
    cat("\nRestricted:", attr(x$x, "is.strat"), "\nFixed margins:", attr(x$x, "fixedmar"))
    if (!is.na(attr(x$x, "shuffle"))) {
        if (attr(x$x, "shuffle")=="ind") cat("\nIndividuals")
        if (attr(x$x, "shuffle")=="samp") cat("\nSamples")
        if (attr(x$x, "shuffle")=="both") cat("\nIndividuals and samples")
    cat(" are shuffled")
    }
    cat("\n\nMatrix dimensions:", nrow(x$x$orig), "rows,", ncol(x$x$orig), "columns")
    cat("\nSum of original matrix:", sum(x$x$orig))
    cat("\nFill of original matrix:", round(sum(x$x$orig>0)/(nrow(x$x$orig)*ncol(x$x$orig)),digits))
    cat("\nNumber of permuted matrices:", n,"\n")
    cat("\nMatrix sums retained:", round(100 * sum(x$sum) / n, digits), "%")
    cat("\nMatrix fill retained:", round(100 * sum(x$fill) / n, digits), "%")
    cat("\nRow sums retained:   ", round(100 * sum(x$rowsums) / n, digits), "%")
    cat("\nColumn sums retained:", round(100 * sum(x$colsums) / n, digits), "%")
    cat("\nRow incidences retained:   ", round(100 * sum(x$browsums) / n, digits), "%")
    cat("\nColumn incidences retained:", round(100 * sum(x$bcolsums) / n, digits), "%")
    if (!is.null(x$strsum))
        cat("\nSums within strata retained:", round(100 * sum(x$strsum) / n, digits), "%")
    cat("\n\nBray-Curtis dissimilarities among original and permuted matrices:\n")
    print(summary(x$bray))
    cat("\nChi-squared for original matrix: ", round(attr(x$chisq, "chisq.orig"), digits), ")\n", sep = "")
    cat("Chi-squared values among expected and permuted matrices:\n")
    print(summary(x$chisq))
invisible(x)
}
