"AIC.radfit" <-
    function (object, k = 2, ...) 
{
    mods <- object$models
    unlist(lapply(mods, AIC, k = k))
}
