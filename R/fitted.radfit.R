`fitted.radfit` <-
    function(object, ...)
{
    out <- sapply(object$models, fitted)
    if (!length(object$y))
        out <- numeric(length(object$models))
    if (length(object$y) <= 1) 
        out <- structure(as.vector(out), dim = c(1, length(object$models)),
                         dimnames = list(names(object$y), names(object$models)))
    out
}

`fitted.radfit.frame` <-
    function(object, ...)
{
    lapply(object, fitted, ...)
}
