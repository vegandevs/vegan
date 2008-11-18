"fisherfit" <-
    function (x, ...) 
{
    Dev.logseries <- function(n.r, p, N) {
        r <- as.numeric(names(n.r))
        x <- N/(N + p)
        logmu <- log(p) + log(x) * r - log(r)
        lhood <- -sum(n.r * (logmu - log(n.r)) + 1) - p * log(1 - 
                                                              x)
        lhood
    }
    tmp <- as.rad(x)
    N <- sum(x)
    tmp <- tmp/N
    p <- 1/sum(tmp^2)
    n.r <- as.fisher(x)
    LSeries <- nlm(Dev.logseries, n.r = n.r, p = p, N = N, 
                   hessian = TRUE, ...)
    LSeries$df.residual <- sum(x > 0) - 1
    LSeries$nuisance <- N/(N + LSeries$estimate)
    LSeries$fisher <- n.r
    class(LSeries) <- "fisherfit"
    LSeries
}
