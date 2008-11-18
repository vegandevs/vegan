"fitted.procrustes" <-
  function(object, truemean = TRUE, ...)
{
  fit <- object$Yrot
  if (truemean)
    fit <- sweep(fit, 2, object$translation, "+")
  fit
}
