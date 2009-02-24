## S3 print method for permat
`print.permat` <-
function(x, digits=3, ...)
{
    cat("Object of class 'permat' with ", attr(x, "times"), " simulations\n", sep="")
    cat("\nMatrix type:", attr(x, "mtype"), "\nPermutation type:", attr(x, "ptype"))
    if (attr(x, "ptype") == "swap") {
        cat("\nMethod: ", attr(x, "method"), sep = "")
        if (attr(x, "method") != "quasiswap") {
            cat(", burnin: ", attr(x, "burnin"), sep = "")
            cat(", thin: ", attr(x, "thin"), sep = "")
        }
    }
    cat("\nRestricted:", attr(x, "is.strat"), "\nFixed margins:", attr(x, "fixedmar"))
    if (!is.na(attr(x, "shuffle"))) {
        if (attr(x, "shuffle")=="ind") cat("\nIndividuals")
        if (attr(x, "shuffle")=="samp") cat("\nSamples")
        if (attr(x, "shuffle")=="both") cat("\nIndividuals and samples")
    cat(" are shuffled")
    }
    cat("\n")
    invisible(x)
#    cat("\n\nMatrix dimensions:", nrow(x$orig), "rows,", ncol(x$orig), "columns")
#    cat("\nSum of original matrix:", sum(x$orig))
#    cat("\nFill of original matrix:", round(sum(x$orig>0)/(nrow(x$orig)*ncol(x$orig)),digits))
#    cat("\nNumber of permuted matrices:", attr(x, "times"),"\n")
}
