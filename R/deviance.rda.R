`deviance.rda` <-
    function(object, ...)
{
    if (is.null(object$CA))
        0
    else
        object$CA$tot.chi * (nrow(object$CA$Xbar) - 1)
}
