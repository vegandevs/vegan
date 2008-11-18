`print.metaMDS` <-
    function (x, ...) 
{
    cat("\nCall:\n")
    cat(deparse(x$call), "\n\n")
    cat("Nonmetric Multidimensional Scaling using isoMDS (MASS package)\n\n")
    cat("Data:    ", x$data, "\n")
    cat("Distance:", x$distance, "\n\n")
    cat("Dimensions:", x$dims, "\n")
    cat("Stress:    ", x$stress, "\n")
    if (x$converged) 
        cat("Two convergent solutions found after", x$tries, 
            "tries\n")
    else cat("No convergent solutions - best solution after", 
             x$tries, "tries\n")
    z <- x$points
    scal <- c(if (attr(z, "centre")) "centring",
              if (attr(z, "pc")) "PC rotation",
              if (attr(z, "halfchange")) "halfchange scaling")
    if (!length(scal)) 
        scal <- "as is"
    cat("Scaling:", paste(scal, collapse = ", "), "\n")
    if (all(is.na(x$species))) {
        cat("Species: scores missing\n")
    } else {
        spattr <- attr(x$species, "shrinkage")
        if (is.null(spattr))
            cat("Species: non-expanded scores ")
        else
            cat("Species: expanded scores ")
        if(attr(x$species, "old.wa"))
            cat("based on untransformed data\n")
        else
            cat("based on", sQuote(x$data), "\n")
        
    }
    cat("\n")
    invisible(x)
}
