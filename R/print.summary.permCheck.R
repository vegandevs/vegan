`print.summary.permCheck` <- function(x, ...)
{
    cat(paste("Number of possible permutations:", x$n, "\n"))
    print(x$control)
    invisible(x)
}
