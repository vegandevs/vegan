`extractAIC.cca` <-
    function (fit, scale = 0, k = 2, ...)
{
   n <- nobs(fit)
   edf <- n - df.residual(fit)
   RSS <- deviance(fit)
   dev <- if(scale > 0)
       RSS/scale - n
   else n * log(RSS/n)
   c(edf, dev + k*edf)
}
