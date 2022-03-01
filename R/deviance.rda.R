`deviance.rda` <-
    function(object, ...)
{
    object$CA$tot.chi * (nobs(object) - 1)
}
