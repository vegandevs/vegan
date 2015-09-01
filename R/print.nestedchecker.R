"print.nestedchecker" <-
function(x, ...)
{
    cat("Checkerboard Units    :", format(x$statistic), "\n")
    cat("C-score (species mean):", format(x$C.score), "\n")
    invisible(x)
}

