`as.ts.permat` <-
    function(x, type = "bray", ...)
{
    type <- match.arg(type, c("bray", "chisq"))
    out <- summary(x)[[type]]
    if (!is.ts(out)) {
        seqmethods <- sapply(make.commsim(), function(z) make.commsim(z)$isSeq)
        seqmethods <- names(seqmethods)[seqmethods]
        ## seqmethods <- c("swap", "tswap", "abuswap")
        stop(gettextf("as.ts available only for sequential methods %s",
                      paste(seqmethods, collapse=", ")))
    }
    out
}
