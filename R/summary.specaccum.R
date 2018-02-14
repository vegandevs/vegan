`summary.specaccum` <-
    function(object, ...)
{
    if (is.null(object$perm))
        stop("specific summary available only for method=\"random\"")
    else {
        tmp <- summary(t(object$perm), ...)
        colnames(tmp) <- paste(1:ncol(tmp), "sites")
        tmp
    }
}
