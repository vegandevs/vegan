`fitted.radfit` <-
    function(object, ...)
{
    matrix(sapply(object$models, fitted), ncol=length(object$models))
}

`fitted.radfit.frame` <-
    function(object, ...)
{
    lapply(object, fitted, ...)
}
