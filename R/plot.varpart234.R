`plot.varpart234` <-
    function(x, cutoff = 0, digits = 1, ...)
{
    vals <- x$indfract[, 3]
    is.na(vals) <- vals < cutoff
    if (cutoff >= 0)
        vals <- round(vals, digits+1)
    labs <-  format(vals, digits=digits, nsmall=digits+1)
    labs <- gsub("NA", "", labs, fixed = TRUE)
    showvarparts(x$nsets, labs, ...)
    if (anyNA(vals)) {
        localMtext <- function(..., Xnames, cutoff) mtext(...)
        localMtext(paste("Values <", cutoff," not shown", sep=""), side=1, ...)
    }
    invisible()
}

