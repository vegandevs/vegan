SSlomolino <-
    selfStart(~ Asym/(1 + slope^log(xmid/area)),
              function(mCall, data, LHS)
{
    xy <- sortedXyData(mCall[["area"]], LHS, data)
    ## wild guess for starting values: assume Asym is 1.2*max(S), fit
    ## log(S) ~ log(area), estimate xmid as area of that model giving
    ## S/2 and assume slope is 1/slope of that model.
    .S <- max(xy[,"y"])
    .b <- as.vector(coef(lm(log(xy[,"y"]) ~ log(xy[,"x"]))))
    value <- c(.S*1.2, exp(log(.S/2)/.b[2] - .b[1]), 1/.b[2])
    names(value) <- mCall[c("Asym","xmid", "slope")]
    value
},
c("Asym","xmid","slope"))
