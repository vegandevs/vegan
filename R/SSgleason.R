SSgleason <-
    selfStart(~ k + slope*log(area),
              function(mCall, data, LHS)
{
    ## Gleason is a linear model: starting values are final ones
    xy <- sortedXyData(mCall[["area"]], LHS, data)
    value <- as.vector(coef(lm(xy[,"y"] ~ log(xy[,"x"]))))
    names(value) <- mCall[c("k","slope")]
    value
},
c("k","slope"))
