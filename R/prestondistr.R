"prestondistr" <-
    function (x, truncate = -1,  ...) 
{
    fun <- function(par, x, truncate) {
        up <- dnorm(x, par[1], par[2], log = TRUE)
        dn <- pnorm(truncate, par[1], par[2], lower = FALSE)
        -sum(up - log(dn))
    }
    x <- x[x > 0]
    logx <- log2(x)
    p <- c(mean(logx), sd(logx))
    sol <- optim(p, fun, x = logx, truncate = truncate)
    p <- sol$par
    area <- pnorm(truncate, p[1], p[2], lower = FALSE)
    scale <- length(x)/sqrt(2 * pi)/p[2]/area
    p <- c(p, scale)
    oct <- as.preston(x)
    x <- as.numeric(names(oct))
    fit <- p[3] * exp(-(x - p[1])^2/2/p[2]^2) 
    names(p) <- c("mode", "width", "S0")
    out <- list(freq = oct, fitted = fit, coefficients = p)
    out$method <- "maximized likelihood to log2 abundances"
    class(out) <- "prestonfit"
    out
}
