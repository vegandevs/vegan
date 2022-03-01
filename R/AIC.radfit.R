### these functions are defined _ex machina_ for radline objects which
### inherit from glm. Here we define them for radfit objects where
### object$models is a list of radline objects

`AIC.radfit` <-
    function (object, k = 2, ...) 
{
    sapply(object$models, AIC, k = k, ...)
}

`deviance.radfit` <-
    function(object, ...)
{
    sapply(object$models, deviance, ...)
}

`logLik.radfit` <-
    function(object, ...)
{
    sapply(object$models, logLik, ...)
}

### Define also for radfit.frames which are lists of radfit objects

`AIC.radfit.frame` <-
    function(object, k = 2, ...)
{
    sapply(object, AIC, k = k, ...)
}

`deviance.radfit.frame` <-
    function(object, ...)
{
    sapply(object, deviance, ...)
}

`logLik.radfit.frame` <-
    function(object, ...)
{
    sapply(object, logLik, ...)
}
