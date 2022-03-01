`print.fisherfit` <-
    function (x, ...) 
{
    cat("\nFisher log series model\n")
    cat("No. of species:", sum(x$fisher), "\n")
    cat("Fisher alpha:  ", x$estimate, "\n\n")
    invisible(x)
}
