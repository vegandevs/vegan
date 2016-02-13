"humpfit" <-
    function (mass, spno, family = poisson, start) 
{
    hump <- function(p, mass, spno, ...) {
        x <- ifelse(mass < p[1], mass/p[1], p[1] * p[1]/mass/mass)
        fv <- p[3] * log(1 + p[2] * x/p[3])
        n <- wt <- rep(1, length(x))
        dev <- sum(dev.resids(spno, fv, wt))
        aicfun(spno, n, fv, wt, dev)/2
    }
    fam <- family(link = "identity")
    aicfun <- fam$aic
    dev.resids <- fam$dev.resids
    if (missing(start))
        p <- c(mean(mass), 100, 10)
    else
        p <- start
    fit <- nlm(hump, p = p, mass = mass, spno = spno, hessian = TRUE)
    p <- fit$estimate
    names(p) <- c("hump", "scale", "alpha")
    x <- ifelse(mass < p[1], mass/p[1], p[1] * p[1]/mass/mass)
    fv <- p[3] * log(1 + p[2] * x/p[3])
    res <- dev.resids(spno, fv, rep(1, length(x)))
    dev <- sum(res)
    residuals <- spno - fv
    aic <- fit$minimum * 2 + 6
    rdf <- length(x) - 3
    out <- list(nlm = fit, family = fam, y = spno, x = mass, 
                coefficients = p, fitted.values = fv, aic = aic, rank = 3, 
                df.residual = rdf, deviance = dev, residuals = residuals, 
                prior.weights = rep(1, length(x)))
    class(out) <- c("humpfit", "glm")
    out
}
