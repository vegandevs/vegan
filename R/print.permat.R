## S3 print method for permat
`print.permat` <-
function(x, digits=3, ...)
{
    cat("Object of class 'permat' with ", attr(x, "times"), " simulations\n", sep="")
    cat("\nMatrix type:", attr(x, "mtype"), "\nPermutation type:", attr(x, "ptype"))
    cat("\nMethod: ", attr(x, "method"), sep = "")
    if (attr(x, "ptype") == "swap") {
        if (!is.na(attr(x, "burnin")))
            cat(", burnin: ", attr(x, "burnin"), sep = "")
        if (!is.na(attr(x, "thin")))
            cat(", thin: ", attr(x, "thin"), sep = "")
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
}
