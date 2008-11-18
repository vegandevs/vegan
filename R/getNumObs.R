`getNumObs` <- function(object, ...) UseMethod("getNumObs")

`getNumObs.default` <- function(object, ...)
{
    nrow(scores(object, display = "sites"))
}

`getNumObs.numeric` <- function(object, ...)
{
    length(object)
}

`getNumObs.integer` <- function(object, ...)
{
    getNumObs.numeric(object)
}
