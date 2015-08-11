`print.varpart` <-
    function (x, ...)
{
    cat("\nPartition of", x$inert, "in", x$RDA, "\n\n")
    writeLines(strwrap(pasteCall(x$call)))
    if (x$scale)
        cat("Columns of Y were scaled to unit variance\n")
    if (!is.null(x$transfo))
        cat("Species transformation: ", x$transfo)
    cat("\n")
    cat("Explanatory tables:\n")
    cat(paste(paste(paste("X", seq_along(x$tables), sep=""),":  ",
                    x$tables, sep=""), collapse="\n"), "\n\n")
    print(x$part, ...)
    invisible(x)
}

