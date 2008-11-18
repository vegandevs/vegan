"print.specaccum" <-
    function(x, ...)
{
    cat("Species Accumulation Curve\n")
    cat("Accumulation method:", x$method)
    if (x$method == "random") {
        cat(", with ", ncol(x$perm), " permutations", sep="")
    }
    cat("\n")
    cat("Call:", deparse(x$call), "\n\n")
    mat <- rbind(Sites = x$sites, Richness = x$richness, sd=x$sd)
    colnames(mat) <- rep("", ncol(mat))
    print(mat)
    invisible(x)
}
