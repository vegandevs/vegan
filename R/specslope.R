#' The Slope of Species Accumulation Curve at Given Point
#' 
#' Function evaluates the derivative of the species accumulation curve
#' for accumulation methods built upon analytic accumulation
#' methods. These methods are \code{exact}, \code{rarefaction} and
#' \code{coleman}. These methods can be evaluated at any sample size,
#' including non-integer values. For other methods, you must look at
#' the differences between consecutive steps, using
#' \code{diff(predict(mod))}.
#'
#' @param object \code{specaccum} result object fitted with methods
#' \code{"exact"}, \code{"rarefaction"} or \code{"coleman"}.
#' @param at The sample size (number of sites) at which the slope is
#' evaluated. This need not be an integer.
`specslope` <-
    function(object, at)
{
    if (!inherits(object, "specaccum"))
        stop("'object' should be a speccaccum() result")
    accepted <- c("exact", "rarefaction", "coleman")
    if (!(object$method %in% accepted))
        stop("accumulation method must be one of: ",
             paste(accepted, collapse=", "))
    ## The following functions are completely defined by species
    ## frequencies
    f <- object$freq
    n <- length(object$sites)
    switch(object$method,
           exact = {
               d <- digamma(pmax(n-at+1, 1)) - digamma(pmax(n-f-at+1, 1))
               g <- lgamma(pmax(n-f+1,1)) + lgamma(pmax(n-at+1,1)) -
                   lgamma(pmax(n-f-at+1, 1)) - lgamma(n+1)
               d <- d*exp(g)
               sum(d[is.finite(d)])
           },
           rarefaction = {
               ## fractional number of individuals at 'at'
               rareslope(f, at/n*sum(f))
           },
           coleman = {
               sum((1 - at/n)^f*f/(n - at))
           })
}
