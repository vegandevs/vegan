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
    UseMethod("specslope")
}

`specslope.specaccum` <-
    function(object, at)
{
    accepted <- c("exact", "rarefaction", "coleman")
    if (!(object$method %in% accepted))
        stop(gettextf("accumulation method must be one of: %s",
             paste(accepted, collapse=", ")))
    ## Funcions should accept a vector of 'at', but usually they
    ## don't. I don't care to change this, and therefore we check the
    ## input.
    if (length(at) > 1 && object$method %in% c("exact", "coleman"))
        stop("'at' can only have a single value")
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
               ## fractional number of individuals at 'at', and slope
               ## for adding whole site instead of one individual
               rareslope(f, at/n*sum(f)) * sum(f)/n
           },
           coleman = {
               sum((1 - at/n)^f*f/(n - at))
           })
}

## Analytical derivatives for NLS regression models in fitspecaccum

`specslope.fitspecaccum` <-
    function(object, at)
{
    ## functions for single set of fitted parameters. Parameters are
    ## given as a single vector 'p' as returned by coef(). Below a
    ## table of original names of 'p':

    ## arrhenius, gitay, gleason: k slope
    ## lomolino: Asym xmid slope
    ## asymp: Asym RO lrc
    ## gompertz: Asym b2 b3
    ## michaelis-menten: Vm K (function SSmicmen)
    ## logis: Asym xmid scal
    ## weibull: Asym Drop lrc pwr
    slope <-
        switch(object$SSmodel,
               "arrhenius" = function(x,p) p[1]*x^(p[2]-1)*p[2],
               "gitay" = function(x,p) 2*(p[1]+p[2]*log(x))*p[2]/x,
               "gleason" = function(x,p) p[2]/x,
               "lomolino" = function(x,p) p[1]*p[3]^log(p[2]/x)*log(p[3])/
                   (1+p[3]^log(p[2]/x))^2/x,
               "asymp" = function(x,p) (p[1]-p[2])*exp(p[3]-exp(p[3])*x),
               "gompertz" = function(x,p) -p[1]*p[2]*p[3]^x*
                   log(p[3])*exp(-p[2]*p[3]^x),
               "michaelis-menten" = function(x,p) p[1]*p[2]/(p[2]+x)^2,
               "logis" = function(x,p) p[1]*exp((x-p[2])/p[3])/
                   (1 + exp((x-p[2])/p[3]))^2/p[3],
               "weibull" = function(x, p) p[2]*exp(p[3]-exp(p[3])*x^p[4])*
               x^(p[4]-1)*p[4])
    ## Apply slope with fitted coefficients at 'at'
    p <- coef(object)
    if (is.matrix(p))          # several fitted models
        out <- apply(p, 2, function(i) slope(at, i))
    else                       # single site drops to a vector
        out <- slope(at, p)
    names(out) <- NULL
    out
}
