`fieller.MOStest` <-
    function (object, level = 0.95) 
{
    smodel <- summary(object$mod)
    var <- smodel$cov.scaled
    fam <- family(object$mod)
    od <- smodel$dispersion
    k <- coef(object$mod)
    b2 <- -2 * k[3]
    u <- -k[2]/2/k[3]
    alpha <- (1-level)/2
    limits <- numeric(2)
    names(limits) <- paste(round(100*(c(alpha, 1-alpha)), 1), "%") 
    wvar <- var[2,2] * od
    uvar <- 4 * var[3,3] * od
    vvar <- -2 * var[2,3] * od
    z <- qnorm(1 - alpha)
    g <- z^2 * uvar/b2^2
    if (g >= 1) {
        limits <- c(NA, NA)
    }
    else {
        x <- u - g * vvar/uvar
        f <- z/b2
        s <- sqrt(wvar - 2 * u * vvar + u^2 * uvar - g * (wvar - 
                                                          vvar^2/uvar))
        limits[1] <- (x - f * s)/(1 - g)
        limits[2] <- (x + f * s)/(1 - g)
    }
    limits
}
