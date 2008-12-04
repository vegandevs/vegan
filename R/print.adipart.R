print.adipart <-
function(x, ...)
{
    cat("adipart with", ncol(x$oecosimu$simulated), "simulations\n")
    att <- attributes(x)
    att$names <- att$class <- att$n.levels <- NULL
    cat(" with", paste(names(att), att, collapse=", "))

    cat("\n\n")
    cl <- class(x)
    if (length(cl) > 1 && cl[2] != "list") {
        NextMethod("print", x)
        cat("\n")
    }
    qu <- apply(x$oecosimu$simulated, 1, quantile, probs=c(0.025, 0.5, 0.975))
    m <- cbind("statistic" = x$oecosimu$statistic,
               "z" = x$oecosimu$z, t(qu),
               "Pr(sim.)"=x$oecosimu$pval)
    printCoefmat(m, ...)
    invisible(x)
}
