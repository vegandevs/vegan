`profile.MOStest` <-
    function(fitted, alpha = 0.01, maxsteps = 10, del = zmax/5, ...)
{
    Pnam <- if(fitted$isHump) "hump" else "pit"
    k <- coef(fitted$mod)
    u <- -k[2]/2/k[3]
    n <- length(residuals(fitted$mod))
    std.error <- fieller.MOStest(fitted, level=0.6)
    std.error <- u - std.error[1]
    if (is.na(std.error))
        std.error <- diff(range(model.matrix(fitted$mod)[,2]))
    OrigDev <- deviance(fitted$mod)
    summ <- summary(fitted$mod)
    DispPar <- summ$dispersion
    fam <- family(fitted$mod)
    Y <- fitted$mod$y
    X <- model.matrix(fitted$mod)[,-3]
    Xi <- X
    if (fam$family %in% c("poisson", "binomial", "Negative Binomial")) {
        zmax <- sqrt(qchisq(1 - alpha/2, 1))
        profName <- "z"
    } else {
        zmax <- sqrt(qf(1 - alpha/2, 1, n - 1))
        profName <- "tau"
    }
    zi <- 0
    prof <- vector("list", length=1)
    names(prof) <- Pnam
    uvi <- u
    for (sgn in c(-1, 1)) {
        step <- 0
        z <- 0
        while((step <- step + 1) < maxsteps && abs(z) < zmax) {
            ui <- u + sgn * step * del * std.error
            Xi[,2] <- (X[,2] - ui)^2
            fm <- glm.fit(x = Xi, y = Y, family=fam,
                          control = fitted$mod$control)
            uvi <- c(uvi, ui)
            zz <- (fm$deviance - OrigDev)/DispPar
            z <- sgn * sqrt(zz)
            zi <- c(zi, z)
        }
        si <- order(zi)
        prof[[Pnam]] <- structure(data.frame(zi[si]), names=profName)
        uvi <- as.matrix(uvi)
        colnames(uvi) <- Pnam
        prof[[Pnam]]$par.vals <- uvi[si, , drop=FALSE]
    }
    of <- list()
    of$coefficients <- structure(Pnam, names=Pnam)
    val <- structure(prof, original.fit = of, summary = summ)
    class(val) <- c("profile.MOStest", "profile.glm", "profile")
    val
}

