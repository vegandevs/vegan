"summary.radfit.frame" <-
function (object, ...) 
{
    labels <- names(object)
    for (i in 1:length(labels)) {
        cat("\n***", labels[i], "***\n")
        print(object[[i]], ...)
    }
    invisible(object)
}

