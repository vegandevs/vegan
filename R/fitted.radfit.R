`fitted.radfit` <-
    function(object, ...)
{
    sapply(object$models, fitted)
}

`fitted.radfit.frame` <-
    function(object, ...)
{
    lapply(object, fitted, ...)
}
