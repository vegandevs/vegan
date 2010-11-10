SSlomolino <-
    selfStart(~ Asym/(1 + slope^log(xmid/area)),
              function(mCall, data, LHS)
{
    xy <- sortedXyData(mCall[["area"]], LHS, data)
    ## approximate with Arrhenius model on log-log
    .p <- coef(lm(log(xy[["y"]]) ~ log(xy[["x"]])))
    ## Asym is value at max(x) but > max(y) and xmid is x which gives
    ## Asym/2
    .Smax <- max(xy[["y"]])*1.1
    .S <- exp(.p[1] + log(max(xy[["x"]])) * (.p[2]))
    .S <- max(.S, .Smax)
    .xmid <- exp((log(.S/2) - .p[1])/.p[2])
    ## approximate slope using y = log(Smax/S - 1)
    .y <- log(.S/xy[["y"]] - 1)
    .x <- xy[["x"]]
    .z <- log(.xmid/xy[["x"]])
    .b <- coef(lm(.y ~ .z))
    ## Adjust Asym: half of y = Asym/2 at xmid
    ##.S <- .S * exp(-0.5 * (.b[1]))
    value <- c(.S, .xmid, exp(.b[2]))
    names(value) <- mCall[c("Asym","xmid", "slope")]
    value
},
c("Asym","xmid","slope"))
