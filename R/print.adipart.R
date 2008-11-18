print.adipart <-
function(x, ...)
{
    cat("Object of class 'adipart' for additive diversity partitioning\n\nCall: ")
    print(x$input$call)
    cat("Method:", attr(x, "method"), "\nDesign:", attr(x, "design"), "\nWeights:", attr(x, "weights"))
    cat("\nIndex:", attr(x, "index"))
    cat("\nNumber of levels:", ncol(x$input$f), "\nNumber of permutations:", attr(x, "times"), "\n")
}
