"print.decorana" <-
function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\nCall:\n")
    cat(deparse(x$call), "\n\n")
    cat(ifelse(x$ira, "Orthogonal", "Detrended"), "correspondence analysis")
    cat(ifelse(!x$ira, paste(" with ", x$mk, " segments.\n", 
        sep = ""), ".\n"))
    if (x$iresc) {
        cat("Rescaling of axes with", x$iresc, "iterations")
        if (x$short) 
            cat(", and shortest axis rescaled", x$short)
        cat(".\n")
    }
    if (!is.null(x$v)) 
        cat("Downweighting of rare species from fraction 1/", x$fraction, ".\n", sep="")
    if (!is.null(x$before)) {
        cat("Piecewise transformation of above-zero abundances:\n")
        print(rbind(before = x$before, after = x$after))
    }
    axlen <- NULL
    if (!x$ira && x$iresc) {
        axlen <- apply(x$rproj, 2, max)
    }
    cat("\n")
    print(rbind(Eigenvalues = x$evals, "Decorana values" = x$evals.decorana, 
        "Axis lengths" = axlen), digits = digits)
    cat("\n")
    invisible(x)
}
