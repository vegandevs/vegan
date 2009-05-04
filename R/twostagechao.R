## Chao et al. 2008. Biometrics 64, 1178-1186
twostagechao <- function(x, order=2, N=nrow(x), m=1, nboot=200, subset)
{
    DNAME <- deparse(substitute(x))
    if (N < 2)
        stop("provide 'N' >= 2")
    if (missing(subset))
        subset <- 1:N
    X <- x[subset,]
    if (N != nrow(X) || N != length(subset))
        stop("'N' and 'subset' must conform")
    tc <- colSums(X)
    x <- xx <- decostand(X, "total")
    if (m > 1)
        x <- 1 - (1 - x)^m
    spp <- ncol(x)
    id <- cbind(unlist(lapply(2:spp, function(z) z:spp)),
        rep(1:(spp-1), (spp-1):1))

    FUN <- function(y) {
        a <- sapply(1:spp, function(z) sum(y[,z])^order)
        b <- sapply(1:spp, function(z) sum(y[,z]^order))
        out <- ((1/(N^order - N)) * (sum(a - b))) / ((1/N) * sum(b))
        out
    }

    if (nboot) {
        BOOT <- sapply(1:spp, function(z) rmultinom(nboot, tc[z], xx[,z]))
        BOOT <- array(BOOT, c(N, nboot, spp))
        BOOTC <- sapply(1:nboot, function(z) FUN(BOOT[,z,]))
    } else BOOTC <- NA

    ESTIMATE <- FUN(x)
    SE <- sd(BOOTC) / sqrt(nboot)
    STATISTIC <- c(ESTIMATE, SE)
    names(STATISTIC) <- c("Estimate", "Std. Error")
    PARAMETER <- c(order, N)
    names(PARAMETER) <- c("q", "N")
    SETEXT <- if (nboot)
        paste("(standard error is based on", nboot, "bootstrap samples)")
        else ""
    METHOD <- paste("Chao's two stage similarity index", SETEXT)

    structure(list(statistic = STATISTIC, parameter = PARAMETER, 
        p.value = NULL, method = METHOD, data.name = DNAME, observed = X, 
        expected = NULL, residuals = NULL), class = c("twostagechao", "htest"))
}
