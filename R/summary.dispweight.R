### summary methods extracts dispweight attributes, and prints a table
### of dispersion statistics

`summary.dispweight`  <-
    function(object, ...)
{
    x <- attributes(object)
    class(x) <- "summary.dispweight"
    x
}

`print.summary.dispweight` <-
    function(x, ...)
{
    tab <- with(x, cbind(D, weights, df, p))
    colnames(tab) <- c("Dispersion", "weight", "d.f.", "Pr(Dispersion)")
    printCoefmat(tab, cs.ind = NA)
    if (!is.na(x$nsimul))
        cat(gettextf("Based on %d simulations on %s nullmodel\n",
                     x$nsimul, x$nullmodel))
    invisible(x)
}
