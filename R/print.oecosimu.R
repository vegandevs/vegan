`print.oecosimu` <-
    function(x, ...)
{
    attr(x$oecosimu$method, "permfun") <- NULL
    cat("oecosimu with", ncol(x$oecosimu$simulated), "simulations\n")
    cat("simulation method", x$oecosimu$method)
    if (length(att <- attributes(x$oecosimu$simulated)) > 1) {
        att$dim <- NULL
        cat(" with", paste(names(att), att, collapse=", "))
    }
    alt.char <- switch(x$oecosimu$alternative,
                       two.sided = "not equal to",
                       less = "less than",
                       greater = "greater than")
    cat("\nalternative hypothesis: true mean is", alt.char, "the statistic")
    ## dim attribute is always there, but print all others

    cat("\n\n")
    cl <- class(x)
    if (length(cl) > 1 && cl[2] != "list") {
        NextMethod("print", x)
        cat("\n")
    }
    probs <- switch(x$oecosimu$alternative,
                    two.sided = c(0.025, 0.5, 0.975),
                    less = c(0, 0.5, 0.95),
                    greater = c(0.05, 0.5, 1))
    qu <- apply(x$oecosimu$simulated, 1, quantile, probs=probs)
    m <- cbind("statistic" = x$oecosimu$statistic,
               "z" = x$oecosimu$z, t(qu),
               "Pr(sim.)"=x$oecosimu$pval)
    printCoefmat(m, ...)
    invisible(x)   
}
