"profile.fisherfit" <-
    function (fitted, alpha = 0.01, maxsteps = 20, del = zmax/5, ...) 
{
    Dev.logseries <- function(n.r, p, N) {
        r <- as.numeric(names(n.r))
        x <- N/(N + p)
        logmu <- log(p) + log(x) * r - log(r)
        lhood <- -sum(n.r * (logmu - log(n.r)) + 1) - p * log(1 -
                                                              x)
        lhood
    }
    par <- fitted$estimate
    names(par) <- "alpha"
    std.err <- sqrt(diag(solve(fitted$hessian)))
    minll <- fitted$minimum
    nr <- fitted$fisher
    N <- sum(as.numeric(names(nr)) * nr)
    zmax <- sqrt(qchisq(1 - alpha/2, 1))
    zi <- 0
    bi <- par
    for (sgn in c(-1, 1)) {
        step <- 0
        z <- 0
        b <- 0
        while ((step <- step + 1) < maxsteps && abs(z) < zmax) {
            b <- par + sgn * step * del * std.err
            fm <- Dev.logseries(nr, b, N)
            zz <- 2 * (fm - minll)
            if (zz > -0.001) 
                zz <- max(zz, 0)
            else stop("profiling has found a better solution, so original fit had not converged")
            z <- sgn * sqrt(zz)
            bi <- c(bi, b)
            zi <- c(zi, z)
        }
    }
    si <- order(bi)
    out <- list()
    out$alpha <- data.frame(tau = zi[si], par.vals = bi[si])
    attr(out, "original.fit") <- list(coefficients = par, std.err = std.err)
    class(out) <- c("profile.fisherfit", "profile.glm", "profile")
    out
}
