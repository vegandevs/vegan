fitspecaccum <- function(object, model = "michaelis-menten", method = "random",  ...)
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
    ## Only Michaelis-Menten implemented now: no need to switch() yet.
    mods <- apply(SpeciesRichness, 2, function(y) nls(y ~ SSmicmen(x, Vm, K))) 
    object$fitted <- drop(sapply(mods, fitted))
    object$residuals <- drop(sapply(mods, residuals))
    object$coefficients <- drop(sapply(mods, coef))
    object$models <- mods
    object$call <- match.call()
    class(object) <- c("fitspecaccum", class(object))
    object
}
