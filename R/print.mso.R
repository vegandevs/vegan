`print.mso` <-
    function(x,  digits = max(3, getOption("digits") - 3), ...)
{
    NextMethod("print")
    cat("mso variogram:\n\n")
    print(x$vario, digits = digits, ...)
    if(!is.null(attr(x$vario, "control")))
        cat("\n", howHead(attr(x$vario, "control")), "\n", sep="")
    invisible(x)
}

