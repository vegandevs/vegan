`print.metaMDS` <-
    function (x, ...) 
{
    cat("\nCall:\n")
    cat(deparse(x$call), "\n\n")
    if (x$engine == "monoMDS")
        cat(x$model, "Multidimensional Scaling using monoMDS\n\n")
    else if (x$engine == "isoMDS")
        cat("non-metric Multidimensional Scaling using isoMDS (MASS package)\n\n")
    cat("Data:    ", x$data, "\n")
    cat("Distance:", x$distance, "\n\n")
    cat("Dimensions:", x$ndim, "\n")
    cat("Stress:    ", x$stress, "\n")
    if (inherits(x, "monoMDS")) {
        cat("Stress type", x$isform)
        if(x$model != "linear")
            cat(", ", c("weak", "strong")[x$ities], " ties", sep = "")
        cat("\n")
    }
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
        cat("based on", sQuote(x$data), "\n")
    }
    cat("\n")
    invisible(x)
}
