`residuals.procrustes` <-
  function (object, ...)
{
  distance <- object$X - object$Yrot
  resid <- rowSums(distance^2)
  sqrt(resid)
}
