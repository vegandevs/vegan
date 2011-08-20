fitspecaccum <-
    function(object, model, method = "random",  ...)
{
    MODELS <- c("arrhenius", "gleason", "gitay", "lomolino",
                "asymp", "gompertz", "michaelis-menten", "logis",
                "weibull")
    model <- match.arg(model, MODELS)
    if (!inherits(object, "specaccum")) 
        object <- specaccum(object, method = method, ...)
    if (is.null(object$perm))
        SpeciesRichness <- as.matrix(object$richness)
    else
        SpeciesRichness <- object$perm
    if (!is.null(object$individuals))
        x <- object$individuals
    else
        x <- object$sites
    mods <- switch(model,
        "arrhenius" = apply(SpeciesRichness, 2,
             function(y) nls(y ~ SSarrhenius(x, k, z), ...)),
        "gleason" = apply(SpeciesRichness, 2,
             function(y) nls(y ~ SSgleason(x, k, slope), ...)),
        "gitay" = apply(SpeciesRichness, 2,
             function(y) nls(y ~ SSgitay(x, k, slope), ...)),
        "lomolino" = apply(SpeciesRichness, 2,
             function(y) nls(y ~ SSlomolino(x, Asym, xmid, slope), ...)),
        "asymp" = apply(SpeciesRichness, 2,
             function(y) nls(y ~ SSlogis(x, Asym, xmid, scal), ...)),
        "gompertz" = apply(SpeciesRichness, 2,
             function(y) nls(y ~ SSgompertz(x, Asym, xmid, scal), ...)),
        "michaelis-menten" = apply(SpeciesRichness, 2,
             function(y) nls(y ~ SSmicmen(x, Vm, K), ...)),
        "logis" = apply(SpeciesRichness, 2,
             function(y) nls(y ~ SSlogis(x, Asym, xmid, scal), ...)),
        "weibull" = apply(SpeciesRichness, 2,
             function(y) nls(y ~ SSweibull(x, Asym, Drop, lrc, par), ...))
                   )
    object$fitted <- drop(sapply(mods, fitted))
    object$residuals <- drop(sapply(mods, residuals))
    object$coefficients <- drop(sapply(mods, coef))
    object$models <- mods
    object$call <- match.call()
    class(object) <- c("fitspecaccum", class(object))
    object
}

### plot function

`plot.fitspecaccum` <-
    function(x, col = par("fg"), lty = 1, 
             xlab = "Sites", ylab = x$method, ...)
{
    fv <- fitted(x)
    matplot(x$sites, fv, col = col, lty = lty, pch = NA,
            xlab = xlab, ylab = ylab, type = "l", ...)
    invisible()
}
