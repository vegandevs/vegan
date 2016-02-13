`extractAIC.cca` <-
    function (fit, scale = 0, k = 2, ...)
{
   n <- nobs(fit)
   edf <- 1
   if (!is.null(fit$CCA$rank)) edf <- edf + fit$CCA$qrank
   if (!is.null(fit$pCCA$rank)) edf <- edf + fit$pCCA$rank
   #edf <- n - fit$CA$rank
   RSS <- deviance(fit)
   dev <- if(scale > 0)
       RSS/scale - n
   else n * log(RSS/n)
   c(edf, dev + k*edf)
}
