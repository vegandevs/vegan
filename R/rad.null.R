"rad.null" <-
    function(x,  family=poisson, ...)
{
    fam <- family(link="log")
    aicfun <- fam$aic
    dev.resids <- fam$dev.resids
    x <- as.rad(x)
    nsp <- length(x)
    wt <- rep(1, nsp)
    if (nsp > 0) { 
        fit <- rev(cumsum(1/nsp:1)/nsp) * sum(x)
        aic <- aicfun(x, nsp, fit, wt, deviance)
    }
    else {
        fit <- NA
        aic <- NA
    }
    res <- dev.resids(x, fit, wt)
    deviance <- sum(res)
    residuals <- x - fit
    rdf <- nsp
    p <- NA
    names(p) <- "S"
    out <- list(model = "Brokenstick", family=fam, y = x, coefficients = p,
                fitted.values = fit, aic = aic, rank = 0, df.residual = rdf,
                deviance = deviance, residuals = residuals, prior.weights=wt)
    class(out) <- c("radline", "glm")
    out
}
