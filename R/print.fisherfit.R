"print.fisherfit" <-
    function (x, ...) 
{
    cat("\nFisher log series model\n")
    cat("No. of species:", sum(x$fisher), "\n\n")
    out <- cbind(x$estimate, sqrt(diag(solve(x$hessian))))
    colnames(out) <- c("Estimate", "Std. Error")
    rownames(out) <- "alpha"
    printCoefmat(out)
    cat("\n")
    invisible(x)
}
