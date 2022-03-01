"summary.radfit.frame" <-
function (object, ...)
{
    labels <- names(object)
    for (i in seq_along(labels)) {
        cat("\n***", labels[i], "***\n")
        print(object[[i]], ...)
    }
    invisible(object)
}

