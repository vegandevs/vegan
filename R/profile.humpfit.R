"profile.humpfit" <-
    function(fitted, parm=1:3, alpha=0.01, maxsteps = 20, del = zmax/5, ...)
{
    INSERT3 <- function(vec, fix, val) {
        switch(fix,
               c(val, vec),
               c(vec[1], val, vec[2]),
               c(vec, val)
               )
    }
    HUMP <- function(p, mass, spno, fix, val, ...) {
        b <- INSERT3(p, fix, val)
        x <- ifelse(mass < b[1], mass/b[1], b[1]*b[1]/mass/mass)
        fv <- b[3] * log(1 + b[2]*x/b[3])
        n <- wt <- rep(1, length(x))
        dev <- sum(dev.resids(spno, fv, wt))
        aicfun(spno, n, fv, wt, dev)/2
    }
    dev.resids <- fitted$family$dev.resids
    aicfun <- fitted$family$aic
    minll <- fitted$nlm$minimum
    p <- coefficients(fitted)
    pv0 <- t(as.matrix(p))
    n <- length(fitted$residuals)
    Pnames <- names(p)
    summ <- summary(fitted)
    dispersion <- summ$dispersion
    std.err <- summ$est[, "Std. Error"]
    if (summ$family == "poisson") {
        zmax <- sqrt(qchisq(1 - alpha/2, 3))
        profName <- "z"
    } else {
        zmax <- sqrt(3 * qf(1 - alpha/2, 3, n-3))
        profName <- "tau"
    }
    prof <- vector("list", length = length(parm))
    names(prof) <- Pnames[parm]
    for (i in parm) {
        zi <- 0
        par <- pv0
        pvi <- pv0[-i]
        pi <- Pnames[i]
        for (sgn in c(-1, 1)) {
            step <- 0
            z <- 0 
            while ((step <- step + 1) < maxsteps && abs(z) < zmax) {
                bi <- p[i] + sgn * step * del * std.err[i]
                fm <- nlm(HUMP, p = pvi, mass = fitted$x, spno = fitted$y, fix = i, val = bi)
                pvi <- fm$estimate
                ri <- INSERT3(pvi, i, bi)
                names(ri) <- Pnames
                par <- rbind(par, ri)
                zz <- 2*(fm$minimum - minll)/dispersion
                if (zz > -0.001)
                    zz <- max(0, zz)
                else
                    stop("profiling has found a better solution, so original fit had not converged:\n", Pnames[i], ": ", bi)
                z <- sgn*sqrt(zz)
                zi <- c(zi, z)
            }
        }
        si <- order(zi)
        prof[[pi]] <- structure(data.frame(zi[si]), names= profName)
        prof[[pi]]$par.vals <- par[si,]
    }
    val <- structure(prof, original.fit = fitted, summary = summ)
    class(val) <- c("profile.humpfit", "profile.glm", "profile")
    val
}
