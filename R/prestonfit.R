"prestonfit" <-
    function (x, ...) 
{
    x <- as.preston(x)
    oct <- as.numeric(names(x))
    fit <- glm(x ~ oct + I(oct^2), family = poisson)
    fv <- fitted(fit)
    p <- coef(fit)
    if (!is.na(p[3]) && p[3] < 0) {
        mu <- -p[2]/2/p[3]
        sd <- sqrt(-1/2/p[3])
        S0 <- exp(p[1] - p[2]^2/4/p[3])
        p <- c(mu, sd, S0)
    }
    else {
        p <- rep(NA, 3)
    }
    names(p) <- c("mode", "width", "S0")
    out <- list(freq = unclass(x), fitted = fv, coefficients = p)
    out$method = "Poisson fit to octaves"
    class(out) <- c("prestonfit")
    out
}
