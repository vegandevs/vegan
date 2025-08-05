### Method functions for monoMDS objects

`print.monoMDS` <-
    function(x, ...)
{
    cat("\nCall:\n")
    cat(deparse(x$call), "\n\n")
    modlab <- switch(x$model,
                     global = "Non-metric",
                     local = "Local non-metric",
                     linear = "Linear",
                     hybrid = "Hybrid",
                     x$model)
    cat(paste(modlab, "Multidimensional Scaling\n\n"))
    cat(x$nobj, "points")
    cat(", dissimilarity", sQuote(x$distmethod))
    if (!is.null(x$distcall))
        cat(", call", sQuote(x$distcall))
    cat("\n\n")
    cat("Dimensions:", x$ndim, "\n")
    cat("Stress:    ", x$stress, "\n")
    cat("Stress type", x$isform)
    if (x$model != "linear")
        cat(", ", c("weak", "strong")[x$ities], " ties", sep = "")
    cat("\n")
    cat("Scores ")
    if (x$iscal == 1)
        cat("scaled to unit root mean square")
    else
        cat("unscaled")
    if (attr(x$points, "pc"))
        cat(", rotated to principal components")
    cat("\n")
    stoplab <- switch(x$icause,
                      "Maximum number of iterations (maxit) reached",
                      "Stress nearly zero (< smin)",
                      "Stress nearly unchanged (ratio > sratmax)",
                      "Scale factor of gradient nearly zero (< sfgrmin)")
    cat("Stopped after ", x$iters, " iterations: ", stoplab, "\n", sep="")
    invisible(x)
}

`scores.monoMDS` <-
    function(x, display = "sites",  shrink = FALSE, choices, tidy = FALSE, ...)
{
    scores.metaMDS(x, display = display, shrink = shrink, choices, tidy = tidy,
                   ...)
}

`plot.monoMDS` <-
    function(x, display = "sites", choices = c(1,2), type = "t",  ...)
{
    ordiplot(x, display = display, choices = choices, type = type, ...)
}

`points.monoMDS` <-
    function(x, display = "sites", choices = c(1,2), select, ...)
{
    x <- scores(x, display = display, choices = choices)
    if (!missing(select))
        x <- .checkSelect(select, x)
    points (x, ...)
    invisible()
}

`text.monoMDS` <-
    function(x, display = "sites", labels, choices = c(1,2), select, ...)
{
    x <- scores(x, display = display, choices = choices)
    if (!missing(select))
        x <- .checkSelect(select, x)
    if (!missing(labels))
        rownames(x) <- labels
    text.ordiplot(x, what = "sites", labels = rownames(x), ...)
    invisible()
}

