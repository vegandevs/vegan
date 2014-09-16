`print.CCorA` <-
    function(x, ...)
{
    cat("\nCanonical Correlation Analysis\n")
    cat("\nCall:\n")
    cat(deparse(x$call), "\n\n")
    out <- structure(rbind(x$Mat.ranks), dimnames = list("Matrix Ranks", c("Y", "X")))
    print(out, ...)
    cat("\n")
    cat("Pillai's trace: ", format(x$Pillai, ...), "\n")
    cat("\n")
    cat("Significance of Pillai's trace:\n")

    cat("from F-distribution:  ", format.pval(x$p.Pillai), "\n")
    if (x$nperm > 0) {
        cat("based on permutations: ")
        cat(x$p.perm,"\n")
        cat(howHead(x$control), "\n")
    }
    out <- rbind("Eigenvalues" = x$EigenValues, "Canonical Correlations" = x$CanCorr)
    colnames(out) <- colnames(x$Cy)
    printCoefmat(out, ...)
    cat("\n")
    out <- rbind("RDA R squares" = x$RDA.Rsquares, "adj. RDA R squares" = x$RDA.adj.Rsq)
    colnames(out) <- c("Y | X", "X | Y")
    printCoefmat(out, ...)
    cat("\n")
    invisible(x) 
}
