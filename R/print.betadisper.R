`print.betadisper` <- function(x, digits = max(3, getOption("digits") - 3),
                               neigen = 8, ...)
{
    ## limit number of eignvals to neigen
    eig <- eigenvals(x)
    nev <- length(eig)
    ax.lim <- min(nev, neigen)
    ##
    cat("\n")
    writeLines(strwrap("Homogeneity of multivariate dispersions\n",
                       prefix = "\t"))
    cat("\n")
    writeLines(strwrap(pasteCall(x$call)))
    cat(paste("\nNo. of Positive Eigenvalues:", sum(eig > 0)))
    cat(paste("\nNo. of Negative Eigenvalues:", sum(eig < 0)))
    cat("\n\n")
    type <- ifelse(attr(x, "type") == "median", "median", "centroid")
    writeLines(strwrap(paste0("Average distance to ", type, ":\n")))
    print.default(format(x$group.distances, digits = digits), quote = FALSE)
    cat("\n")
    writeLines(strwrap("Eigenvalues for PCoA axes:"))
    if (nev > neigen) {
        writeLines(strwrap(paste0("(Showing ", neigen, " of ", nev,
                                  " eigenvalues)")))
    }
    print.default(format(eig[seq_len(ax.lim)], digits = digits), quote = FALSE)
    invisible(x)
}
