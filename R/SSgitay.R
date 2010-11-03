SSgitay <-
    selfStart(~ (k + slope*log(area))^2,
              function(mCall, data, LHS)
{
    xy <- sortedXyData(mCall[["area"]], LHS, data)
    value <- as.vector(coef(lm(sqrt(xy[,"y"]) ~ log(xy[,"x"]))))
    names(value) <- mCall[c("k","slope")]
    value
},
c("k","slope"))
