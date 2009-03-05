`as.ts.oecosimu` <-
    function(x, ...)
{
    seqmethods <- c("swap", "tswap")
    if (!(x$oecosimu$method %in% seqmethods))
        stop("as.ts available only for sequential methods ",
             paste(seqmethods, collapse=", "))
    as.ts(t(x$oecosimu$simulated))
}
