`print.nestednodf` <-
    function(x, ...)
{    
    cat("N columns  :", format(x$statistic["N.columns"], ...), "\n")
    cat("N rows     :", format(x$statistic["N.rows"], ...), "\n")
    cat("NODF       :", format(x$statistic["NODF"], ...), "\n")
    cat("Matrix fill:", format(x$fill, ...), "\n")
    invisible(x)
}
