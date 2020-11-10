SSarrhenius <-
    selfStart(~ k*area^z,
              function(mCall, data, LHS, ...)
{
    xy <- sortedXyData(mCall[["area"]], LHS, data)
    value <- as.vector(coef(lm(log(pmax(xy[,"y"],1)) ~ log(xy[,"x"]))))
    value[1] <- exp(value[1])
    names(value) <- mCall[c("k","z")]
    value
},
c("k","z"))
