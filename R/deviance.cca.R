`deviance.cca` <-
    function(object, ...)
{
    if (is.null(object$CA))
        0
    else
        object$CA$tot.chi * object$grand.tot
}
