`print.specaccum` <-
    function(x, ...)
{
    cat("Species Accumulation Curve\n")
    cat("Accumulation method:", x$method)
    if (x$method == "random") {
        cat(", with ", ncol(x$perm), " permutations", sep="")
    }
    if (!is.null(x$weights))
        cat(", weighted")
    cat("\n")
    cat("Call:", deparse(x$call), "\n\n")
    mat <- rbind(Sites = x$sites, Individuals = x$individuals, Effort = x$effort,
                 Richness = x$richness, sd=x$sd)
    colnames(mat) <- rep("", ncol(mat))
    print(zapsmall(mat))
    invisible(x)
}
