`fitted.radfit` <-
    function(object, ...)
{
    matrix(sapply(object$models, fitted), ncol=length(object$models))
}
