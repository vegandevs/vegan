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
    hasWeights <- !is.null(object$weights)
    ## scale weights to correspond to the no. of sites
    if (hasWeights) {
        w <- as.matrix(object$weights)
        n <- nrow(w)
        w <- sweep(w, 2, w[n,], "/") * n
    }
    NLSFUN <- function(y, x, model, ...) {
        switch(model,
        "arrhenius" = nls(y ~ SSarrhenius(x, k, z),  ...),
        "gleason" = nls(y ~ SSgleason(x, k, slope),  ...),
        "gitay" = nls(y ~ SSgitay(x, k, slope), ...),
        "lomolino" = nls(y ~ SSlomolino(x, Asym, xmid, slope), ...),
        "asymp" = nls(y ~ SSasymp(x, Asym, R0, lrc), ...),
        "gompertz" = nls(y ~ SSgompertz(x, Asym, xmid, scal), ...),
        "michaelis-menten" = nls(y ~ SSmicmen(x, Vm, K),  ...),
        "logis" = nls(y ~ SSlogis(x, Asym, xmid, scal),  ...),
        "weibull" = nls(y ~ SSweibull(x, Asym, Drop, lrc, par), ...))
    }
    mods <- lapply(seq_len(NCOL(SpeciesRichness)),
                  function(i, ...)
                   NLSFUN(SpeciesRichness[,i],
                          if (hasWeights) w[,i] else x,
                          model, ...), ...)
    object$fitted <- drop(sapply(mods, fitted))
    object$residuals <- drop(sapply(mods, residuals))
    object$coefficients <- drop(sapply(mods, coef))
    object$SSmodel <- model
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
    if (is.null(x$weights))
        fv <- fitted(x)
    else
        fv <- sapply(x$models, predict, newdata = list(x = x$sites))
    matplot(x$sites, fv, col = col, lty = lty, pch = NA,
            xlab = xlab, ylab = ylab, type = "l", ...)
    invisible()
}

`lines.fitspecaccum` <-
    function(x, col = par("fg"), lty = 1, ...)
{
    if (is.null(x$weights))
        fv <- fitted(x)
    else
        fv <- sapply(x$models, predict, newdata= list(x = x$sites))
    matlines(x$sites, fv, col = col, lty = lty, pch = NA, type = "l", ...)
    invisible()
}
