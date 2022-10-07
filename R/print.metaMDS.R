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
    if (x$converged) {
        cat(sprintf(ngettext(x$converged,
                     "Best solution was repeated %d time in %d tries\n",
                     "Best solution was repeated %d times in %d tries\n"),
                     x$converged, x$tries))
    } else {
        cat("Best solution was not repeated after", x$tries, "tries\n")
    }
    cat("The best solution was from try", x$bestry)
    if (x$bestry == 0)
        cat(" (metric scaling or null solution)\n")
    else
        cat(" (random start)\n")
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
        spdata <- attr(x$species, "data")
        if (is.null(spdata))
            spdata <- x$data
        if (is.null(spattr))
            cat("Species: non-expanded scores ")
        else
            cat("Species: expanded scores ")
        cat("based on", sQuote(spdata), "\n")
    }
    cat("\n")
    invisible(x)
}
