`print.oecosimu` <-
    function(x, ...)
{
    xx <- x ## return unmodified input object
    attr(x$oecosimu$method, "permfun") <- NULL
    cat(as.character(attr(x,"call")[[1]]), "object\n\n")
    writeLines(strwrap(pasteCall(attr(x, "call"))))
    cat("\n")
    cat("nullmodel method", sQuote(x$oecosimu$method), "with", 
        ncol(x$oecosimu$simulated), "simulations\n")
    if (length(att <- attributes(x$oecosimu$simulated)) > 1) {
        att$dim <- NULL
        cat("options: ", paste(names(att), att, collapse=", "))
    }
    alt.char <- switch(x$oecosimu$alternative,
                       two.sided = "less or greater than",
                       less = "less than",
                       greater = "greater than")
    cat("\nalternative hypothesis: statistic is", alt.char, "simulated values")
    ## dim attribute is always there, but print all others

    cat("\n\n")

    if (!inherits(x, c("adipart", "hiersimu", "multipart")) &&
        !inherits(x$statistic, c("numeric", "list"))) {
            print(x$statistic)
            cat("\n")
    }
    probs <- switch(x$oecosimu$alternative,
                    two.sided = c(0.025, 0.5, 0.975),
                    greater = c(0.5, 0.95),
                    less = c(0.05, 0.5))
    qu <- apply(x$oecosimu$simulated, 1, quantile, probs=probs, na.rm = TRUE)
    m <- cbind("statistic" = x$oecosimu$statistic,
               "z" = x$oecosimu$z, "mean" = x$oecosimu$means, t(qu),
               "Pr(sim.)"=x$oecosimu$pval)
    printCoefmat(m, cs.ind = 3:(ncol(m)-1), ...)
    if (any(is.na(x$oecosimu$simulated))) {
        nacount <- rowSums(is.na(x$oecosimu$simulated))
        cat("\nNumber of NA cases removed from simulations:\n",
            nacount, "\n")
    }
    invisible(xx)   
}


