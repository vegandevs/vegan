`print.betadisper` <- function(x, digits = max(3, getOption("digits") - 3),
                               ...)
{
    ## limit number of eignvals to 8
    ax.lim <- 8
    ##
    cat("\n")
    writeLines(strwrap("Homogeneity of multivariate dispersions\n",
                       prefix = "\t"))
    cat("\n")
    writeLines(strwrap(pasteCall(x$call)))
    cat(paste("\nNo. of Positive Eigenvalues:", sum(x$eig > 0)))
    cat(paste("\nNo. of Negative Eigenvalues:", sum(x$eig < 0)))
    cat("\n\n")
    type <- ifelse(isTRUE(all.equal(attr(x, "type"), "median")),
                   "medoid", "centroid")
    writeLines(strwrap(paste0("Average distance to ", type, ":\n")))
    print.default(tapply(x$distances, x$group, mean), digits = digits)
    cat("\n")
    writeLines(strwrap("Eigenvalues for PCoA axes:\n"))
    print.default(round(x$eig[seq_len(ax.lim)], digits = digits))
    invisible(x)
}
