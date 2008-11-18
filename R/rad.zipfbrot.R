"rad.zipfbrot" <-
    function (x, family = poisson, ...) 
{
    mandelfun <- function(p, x, ...) {
        brnk <- log(rnk + exp(p))
        sol <- glm(x ~ brnk + offset(off), family = family(link = "log"))
        -logLik(sol)
    }
    x <- as.rad(x)
    rnk <- seq(along = x)
    off <- rep(log(sum(x)), length(x))
    p <- 0
    fam <- family(link = "log")
    if (length(x) > 2) 
        nl <- try(nlm(mandelfun, p = p, x = x, rnk = rnk, off = off, 
                      family = fam, hessian = TRUE, ...))
    if (length(x) < 3) {
        aic <- NA
        dev <- rdf <-  0
        ln <- nl <- NA
        p <- rep(NA, 3)
        fit <- x
        res <- rep(0, length(x))
        wts <- rep(1, length(x))
    }
    else if (inherits(nl, "try-error")) {
        aic <- rdf <- ln <- nl <- dev <-  NA
        p <- rep(NA, 3)
        fit <- res <- wts <- rep(NA, length(x))
    }
    else {
        ln <- glm(x ~ log(rnk + exp(nl$estimate)) + offset(off), 
                  family = family(link = "log"))
        fit <- fitted(ln)
        p <- c(coef(ln), exp(nl$estimate))
        p[1] <- exp(p[1])
        aic <- AIC(ln) + 2
        rdf <- df.residual(ln) - 1
        dev <- deviance(ln)
        res <- ln$residuals
        wts <- weights(ln)
    }
    names(p) <- c("c", "gamma", "beta")
    out <- list(model = "Zipf-Mandelbrot", family = fam, 
                y = x, coefficients = p, fitted.values = fit, aic = aic, 
                rank = 3, df.residual = rdf, deviance = dev, 
                residuals = res, prior.weights = wts)
    class(out) <- c("radline", "glm")
    out
}
