`print.mso` <-
    function(x,  digits = max(3, getOption("digits") - 3), ...)
{
    NextMethod("print", x, digits = digits, ...)
    cat("mso variogram:\n\n")
    print(x$vario, digits = digits, ...)
    invisible(x)
}

