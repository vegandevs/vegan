`as.ts.oecosimu` <-
    function(x, ...)
{
    if  (!x$oecosimu$isSeq)
        stop("as.ts available only for sequential methods")
    startval <- attr(x$oecosimu$simulated, "burnin") + 1 
    thin <- attr(x$oecosimu$simulated, "thin")
    out <- ts(t(x$oecosimu$simulated), start = startval, deltat=thin,
        names = names(x$oecosimu$z))
    attr(out, "burnin") <- NULL
    attr(out, "thin") <- NULL
    out
}
