### na.action = na.exclude puts NA values to removed observations, but
### in constrained ordination WA scores can be found for observations
### with NA values in constraints.
`ordiNAexclude` <-
    function(object, newdata)
{
    ## Embed NA for excluded cases
    nas <- object$na.action
    object$rowsum <- napredict(nas, object$rowsum)
    object$CCA$u <- napredict(nas, object$CCA$u)
    object$CCA$u.eig <- napredict(nas, object$CCA$u.eig)
    object$CCA$wa <- napredict(nas, object$CCA$wa)
    object$CCA$wa.eig <- napredict(nas, object$CCA$wa.eig)
    object$CA$u <- napredict(nas, object$CA$u)
    object$CA$u.eig <- napredict(nas, object$CA$u.eig)
    ## Estimate WA scores for NA cases with newdata of excluded
    ## obseravations
    wa <- predict(object, newdata = newdata, type = "wa", model = "CCA")
    wa.eig <- sweep(wa, 2, sqrt(object$CCA$eig), "*")
    object$CCA$wa[nas,] <- wa
    object$CCA$wa.eig[nas,] <- wa.eig
    wa <- predict(object, newdata = newdata, type = "wa", model = "CA")
    wa.eig <- sweep(wa, 2, sqrt(object$CA$eig), "*")
    object$CA$u[nas,] <- wa
    object$CA$u.eig[nas,] <- wa.eig
    object
}
