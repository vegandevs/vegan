fitspecaccum <-
    function(object, model = c("michaelis-menten", "arrhenius"),
             method = "random",  ...)
{
    model <- match.arg(model)
    if (!inherits(object, "specaccum")) 
        object <- specaccum(object, method = method, ...)
    if (is.null(object$perm))
        SpeciesRichness <- as.matrix(object$richness)
    else
        SpeciesRichness <- object$perm
    if (!is.null(object$inviduals))
        x <- object$individuals
    else
        x <- object$sites
    mods <- switch(model,
        "michaelis-menten" = apply(SpeciesRichness, 2,
             function(y) nls(y ~ SSmicmen(x, Vm, K))),
        "arrhenius" = apply(SpeciesRichness, 2,
                   function(y) nls(y ~ SSarrhenius(x, k, z))))
    object$fitted <- drop(sapply(mods, fitted))
    object$residuals <- drop(sapply(mods, residuals))
    object$coefficients <- drop(sapply(mods, coef))
    object$models <- mods
    object$call <- match.call()
    class(object) <- c("fitspecaccum", class(object))
    object
}
