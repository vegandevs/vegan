`summary.poolaccum` <-
    function(object, display, alpha = 0.05, ...)
{
    probs <- c(alpha/2, 1-alpha/2)
    if (inherits(object, "estaccumR"))
        dislabels <- c("S", "chao", "ace")
    else
        dislabels <- c("S", "chao", "jack1", "jack2", "boot")
    disnames <- colnames(object$means[,-1])
    names(disnames) <- dislabels
    if (missing(display))
        display <- dislabels
    else
        display <- match.arg(display, dislabels, several.ok = TRUE)
    out <- list()
    for (item in display) {
        out[[item]] <- cbind(`N` = object$N,
                             `Mean` = object$means[,disnames[item], drop=FALSE],
                             t(apply(object[[item]], 1, quantile, probs=probs)),
                             `Std.Dev` = apply(object[[item]], 1, sd))
    }
    class(out) <- "summary.poolaccum"
    out
}
