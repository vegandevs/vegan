`print.summary.cca` <-
    function (x, digits = x$digits, ...)
{
    cat("\nCall:\n")
    cat(deparse(x$call), "\n")
    statnam <- if (x$method == "cca")
        "averages"
    else "sums"
    cat("\nPartitioning of ", x$inertia, ":\n", sep = "")
    out <- c(Total = x$tot.chi, Conditioned = x$partial.chi,
             Constrained = x$constr.chi, Unconstrained = x$unconst.chi)
    out <- cbind(Inertia = out, Proportion = out/out[1])
    print(out, digits = digits, ...)
    cat("\nEigenvalues, and their contribution to the", x$inertia,
        "\n")
    if (!is.null(x$partial.chi)) {
        cat("after removing the contribution of conditiniong variables\n")
    }
    cat("\n")
    print(x$cont$importance, ...)
    if (!is.null(x$concont)) {
        cat("\nAccumulated constrained eigenvalues\n")
        print(x$concont$importance, ...)
    }
    cat("\n")
    invisible(x)
}
