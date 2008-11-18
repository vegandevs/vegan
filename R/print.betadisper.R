`print.betadisper` <- function(x, digits = max(3, getOption("digits") - 3),
                             ...)
{
    cat("\n")
    writeLines(strwrap("Homogeneity of multivariate dispersions\n",
                       prefix = "\t"))
    cat("\n")
    writeLines(strwrap(pasteCall(x$call)))
    cat(paste("\nNo. of Positive Eigenvalues:", sum(x$eig > 0)))
    cat(paste("\nNo. of Negative Eigenvalues:", sum(x$eig < 0)))
    cat("\n\n")
    writeLines(strwrap("Average distance to centroid:\n"))
    print.default(tapply(x$distances, x$group, mean), digits = digits)
    cat("\n")
    writeLines(strwrap("Eigenvalues for PCoA axes:\n"))
    print.default(round(x$eig, digits = digits))
    invisible(x)
}
