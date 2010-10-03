SSarrhenius <-
    selfStart(function(area, k, z) k*area^z,
              function(mCall, data, LHS)
{
    xy <- sortedXyData(mCall[["area"]], LHS, data)
    value <- as.vector(coef(lm(log(xy[,"y"]) ~ log(xy[,"x"]))))
    names(value) <- mCall[c("k","z")]
    value
},
c("k","z"))
