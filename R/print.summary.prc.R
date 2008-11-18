"print.summary.prc" <-
function(x, ...)
{
    cat("\nCall:\n")
    cat(deparse(x$call), "\n")
    cat("Species scores:\n")
    print(x$sp, digits=x$digits, ...)
    cat("\nCoefficients for", paste(x$names, collapse=":"), "interaction\n")
    cat(paste("which are contrasts to", x$names[2], x$corner, "\n"))
    cat(paste(c("rows are",", columns are"), x$names[2:1], collapse=""))
    cat("\n")
    print(coef(x), digits = x$digits, ...)
    invisible(x)
}

