`MOStest` <-
    function(x, y, interval, ...)
{
    if (!missing(interval))
        interval <- sort(interval)
    x <- eval(x)
    m0 <- glm(y ~ x + I(x^2), ...)
    k <- coef(m0)
    isHump <- unname(k[3] < 0)
    hn <- if(isHump) "hump" else "pit"
    hump <- unname(-k[2]/2/k[3])
    if (missing(interval))
        p1 <- min(x)
    else
        p1 <- interval[1]
    if (missing(interval))
        p2 <- max(x)
    else
        p2 <- interval[2]
    test <- if (m0$family$family %in% c("binomial", "poisson")) "Chisq" else "F"
    tmp <- glm(y ~ I(x^2 - 2*x*p1), ...)
    ## Chisq test has one column less than F test: extract statistic
    ## and its P value
    statmin <- anova(tmp, m0, test = test)[2, (5:6) - (test == "Chisq")]
    tmp <- glm(y ~ I(x^2 - 2*x*p2), ...)
    statmax <- anova(tmp, m0, test = test)[2, (5:6) - (test == "Chisq")]
    comb <- 1 - (1-statmin[2])*(1-statmax[2])
    comb <- unlist(comb)
    stats <- rbind(statmin, statmax)
    rownames(stats) <- paste(hn, c("at min", "at max"))
    stats <- cbind("min/max" = c(p1,p2), stats)
    stats <- rbind(stats, "Combined" = c(NA, NA, comb))
    vec <- c(p1, p2, hump)
    names(vec) <- c("min", "max", hn)
    vec <- sort(vec)
    isBracketed <- names(vec)[2] == hn
    out <- list(isHump = isHump, isBracketed = isBracketed,
                hump = vec, family = family(m0), coefficients = stats,
                mod = m0)
    class(out) <- "MOStest"
    out
}
