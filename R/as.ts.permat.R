`as.ts.permat` <-
    function(x, type = "bray", ...)
{
    type <- match.arg(type, c("bray", "chisq"))
    out <- summary(x)[[type]]
    if (!is.ts(out)) {
        seqmethods <- c("swap", "tswap", "abuswap")
        stop("as.ts available only for sequential methods ",
            paste(seqmethods, collapse=", "))
    } 
    out
}
