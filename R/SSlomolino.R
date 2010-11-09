SSlomolino <-
    selfStart(~ Asym/(1 + slope^log(xmid/area)),
              function(mCall, data, LHS)
{
    xy <- sortedXyData(mCall[["area"]], LHS, data)
    ## assume Asym is const*max(S)
    .S <- max(xy[["y"]])*1.5
    ## approximate using y = log(Smax/S - 1)
    .y <- log(.S/xy[["y"]] - 1)
    .x <- xy[["x"]]
    ## xmid is the inflection point: fit a parabola to log and find it
    ## there
    .p <- coef(lm(.y ~ .x + I(.x^2)))
    .xmid <- -(.p[2])/2/.p[3] - sqrt(abs(1/2/.p[3]))
     ## estimate slope assuming Asym and xmid are known
    .z <- log(.xmid/xy[["x"]])
    .b <- coef(lm(.y ~ .z))
    ## Adjust Asym: half of y = Asym/2 at xmid
    .S <- .S * exp(-0.5 * (.b[1]))
    value <- c(.S, .xmid, exp(.b[2]))
    names(value) <- mCall[c("Asym","xmid", "slope")]
    value
},
c("Asym","xmid","slope"))
