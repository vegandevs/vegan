"deviance.rda" <-
function(object, ...) object$CA$tot.chi * (nrow(object$CA$Xbar) - 1)
