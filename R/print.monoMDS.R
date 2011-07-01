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
                     hybrid = "Hybrid")
    cat(paste(modlab, "Multidimensional Scaling\n"))
    cat("Dimensions:", x$ndim, "\n")
    cat("Stress:    ", x$stress)
    if (x$isform == 2)
        cat(" (type 2)")
    cat("\n")
    stoplab <- switch(x$icause,
                      "Maximum number of iteration reached",
                      "Stress nearly zero",
                      "Stress nearly unchanged",
                      "Scale factor of gradient nearly zero")
    cat("Stopped after ", x$iters, " iterations: ", stoplab, "\n", sep="")
    invisible(x)
}

`scores.monoMDS` <-
    function(x, ...)
{
    x$points
}

`plot.monoMDS` <-
    function(x, choices = c(1,2), type = "t",  ...)
{
    ordiplot(x, display = "sites", choices = choices, type = type, ...)
}

