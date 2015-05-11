#' Slope of Rarefunction Curve at Given Sample Size
#'
#' Function evaluates the derivative of the rarefaction
#' function at given sample size. The derivative was
#' directly derived from the expression used in \code{rarefy}.
#'
#' @param x Community counts; a single integer vector
#' @param sample Sample size where the derivative is evaluated; can be real
#'
`rareslope` <-
    function(x, sample)
{
    x <- x[x>0]
    J <- sum(x)
    d <- digamma(pmax(J-sample+1, 1)) - digamma(pmax(J-x-sample+1, 1))
    g <- lgamma(pmax(J-x+1, 1)) + lgamma(pmax(J-sample+1, 1)) -
        lgamma(pmax(J-x-sample+1, 1)) - lgamma(J+1)
    d <- d*exp(g)
    sum(d[is.finite(d)])
}
