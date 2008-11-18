"residuals.procrustes" <-
  function (object, ...) 
{
  distance <- object$X - object$Yrot
  resid <- apply(distance^2, 1, sum)
  resid <- sqrt(resid)
  resid
}
