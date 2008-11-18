"print.varpart" <-
function (x, ...) 
{
    cat("\nPartition of variation in RDA\n\n")
    cat("Call:\n")
    cat(deparse(x$call), "\n")
    if (x$scale)
        cat("Columns of Y were scaled to unit variance\n")
    if (!is.null(x$transfo))
        cat("Species transformation: ", x$transfo)
    cat("\n")
    cat("Explanatory tables:\n")
    cat(paste(paste(paste("X", 1:length(x$tables), sep=""),":  ",
                    x$tables, sep=""), collapse="\n"), "\n\n")
    print(x$part, ...)
    invisible(x)
}

