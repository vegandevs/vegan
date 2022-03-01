"metaMDSredist" <-
function(object, ...)
{
    if (!inherits(object, "metaMDS"))
        stop("Needs a metaMDS result object")
    call <- object$call
    call[[1]] <- as.name("metaMDSdist")
    eval(call, parent.frame())
}

