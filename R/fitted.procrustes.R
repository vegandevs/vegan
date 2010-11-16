`fitted.procrustes` <-
  function(object, truemean = TRUE, ...)
{
  fit <- object$Yrot
  if (truemean)
    fit <- sweep(fit, 2, object$xmean, "+")
  fit
}

## Like above, except when takes 'newata'

`predict.procrustes` <-
    function(object, newdata, truemean = TRUE, ...)
{
    if (missing(newdata))
        return(fitted(object, truemean = truemean))
    Y <- as.matrix(newdata)
    ## scaling and rotation
    Y <- object$scale * Y %*% object$rotation
    ## translation
    if (truemean)
        Y <- sweep(Y, 2, object$translation, "+")
    Y
}
