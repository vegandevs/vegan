## S3 print method for summary.permat
`print.summary.permat` <-
function(x, digits=2, ...)
{
    bray <- x$bray
    restr <- x$restr
    test <- x$test
    x <- x$x
    cat("Summary of object of class 'permat'\n\nCall: ")
    print(x$call)
    cat("\nMatrix type:", attr(x, "mtype"), "\nPermutation type:", attr(x, "ptype"))
    if (attr(x, "ptype") == "swap")
        cat("\nMethod:", attr(x, "method"))
    cat("\nRestricted:", restr, "\nFixed margins:", attr(x, "fixedmar"))
    if (!is.na(attr(x, "shuffle"))) {
        if (attr(x, "shuffle")=="ind") cat("\nIndividuals")
        if (attr(x, "shuffle")=="samp") cat("\nSamples")
        if (attr(x, "shuffle")=="both") cat("\nIndividuals and samples")
    cat(" are shuffled")}
    cat("\n\nMatrix dimensions:", nrow(x$orig), "rows,", ncol(x$orig), "columns")
    cat("\nSum of original matrix:", sum(x$orig))
    cat("\nFill of original matrix:", round(sum(x$orig>0)/(nrow(x$orig)*ncol(x$orig)),digits))
    cat("\nNumber of permuted matrices:", attr(x, "times"),"\n")
    cat("\nMatrix sums retained:", round(100*test[1], digits), "%")
    cat("\nMatrix fill retained:", round(100*test[2], digits), "%")
    cat("\nRow sums retained:   ", round(100*test[3], digits), "%")
    cat("\nColumn sums retained:", round(100*test[4], digits), "%")
    if (restr) cat("\nSums within strata retained:", round(100*test[5], digits), "%")
    cat("\n\nBray-Curtis dissimilarities among original and permuted matrices:\n")
    print(summary(bray))
invisible(NULL)
}
