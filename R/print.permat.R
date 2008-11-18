## S3 print method for permat
`print.permat` <-
function(x, digits=3, ...)
{
    if (attr(x, "ptype") != "sar" & !is.null(x$specs$reg) | !is.null(x$specs$hab))
        restr <- TRUE else restr <- FALSE
    cat("Object of class 'permat'\n\nCall: ")
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
}
