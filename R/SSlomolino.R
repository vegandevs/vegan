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
    ## approximate slope for log(Asym/y - 1) ~ log(xmid/x) + 0
    .y <- log(.S/xy[["y"]] - 1)
    .z <- log(.xmid/xy[["x"]])
    .b <- coef(lm(.y ~ .z - 1))
    value <- c(.S, .xmid, exp(.b))
    names(value) <- mCall[c("Asym","xmid", "slope")]
    value
},
c("Asym","xmid","slope"))
