`print.MOStest` <-
    function(x, ...)
{
    cat("\nMitchell-Olds and Shaw test\n")
    cat("Null: hump of a quadratic linear predictor is at min or max\n")
    print(x$family)
    print(x$hump)
    if (!x$isBracketed)
        cat("***** Caution: hump/pit not bracketed by the data ******\n")
    cat("\n")
    printCoefmat(coef(x), has.P=TRUE, na.print="")
    invisible(x)
}
