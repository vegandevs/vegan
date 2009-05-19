print.adipart <-
function(x, ...)
{
    n <- if (is.null(x$oecosimu$simulated))
        0 else ncol(x$oecosimu$simulated)
    cat("adipart with", n, "simulations\n")
    att <- attributes(x)
    att$names <- att$call <- att$class <- att$n.levels <- att$terms <- att$model <- NULL
    cat("with", paste(names(att), att, collapse=", "))

    cat("\n\n")
    cl <- class(x)
    if (length(cl) > 1 && cl[2] != "list") {
        NextMethod("print", x)
        cat("\n")
    }
    if (!is.null(x$oecosimu$simulated)) {
        tmp <- x$oecosimu$simulated
    } else {
        tmp <- data.matrix(x$oecosimu$statistic)
    }
    qu <- apply(tmp, 1, quantile, probs=c(0.025, 0.5, 0.975))
    m <- cbind("statistic" = x$oecosimu$statistic,
               "z" = x$oecosimu$z, t(qu),
               "Pr(sim.)"=x$oecosimu$pval)
    printCoefmat(m, ...)
    invisible(x)
}
